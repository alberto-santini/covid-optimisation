import pandas as pd
import geopandas as gpd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from os import system
from datetime import datetime, timedelta
from math import radians, cos, sin, asin, sqrt

###########################################################################
# Utilities
import json

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)

def haversine(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon/2) ** 2
    return 6371 * 2 * asin(sqrt(a))
#
###########################################################################

###########################################################################
# Instance destination path:
inst_path = '/home/alberto/local/src/my-github/covid/data/swab-test/generated'
#
###########################################################################

###########################################################################
# Utility to read real data from Italy's DPC
def get_regional_data(start_date='2020-04-01'):
    """ Gets regional data from the dataset of Italian's Civil Defence and
        adds geographic coordinates, from a shapefile with Italian region.

        Parameters:
            - start_date: First day to consider, format: 'yyyy-mm-dd'

        Returns a pair of:
            - a dictionary, whose keys are region names and values are data
              relative to that region
            - a pandas dataframe with data for all regions
    """

    italy = gpd.read_file('gadm36_ITA_1.shp')
    start_date = datetime.strptime(start_date, "%Y-%m-%d").date()
    
    regional_data = dict()
    r = pd.read_csv('dpc-covid19-ita-regioni.csv', parse_dates=['data'])
    r['data'] = pd.to_datetime(r.data).dt.date
    
    for reg in r.denominazione_regione.unique():
        rr = r[r.denominazione_regione == reg].copy()
        rr['nuovi_tamponi'] = rr.tamponi.diff(periods=1).clip(lower=0)
        rr = rr[rr.data >= start_date]
        rr = rr[['data', 'denominazione_regione', 'nuovi_tamponi', 'nuovi_positivi']].copy()
    
        regional_data[reg] = rr
        
    # Merge Trento and Bolzano
    regional_data['Trentino Alto Adige'] = pd \
        .concat([regional_data['P.A. Bolzano'], regional_data['P.A. Trento']]) \
        .groupby('data') \
        .sum() \
        .reset_index()
    regional_data.pop('P.A. Bolzano')
    regional_data.pop('P.A. Trento')
    
    # Rename regions
    regional_data['Emilia Romagna'] = regional_data.pop('Emilia-Romagna')
    regional_data['Valle Aosta'] = regional_data.pop('Valle d\'Aosta')
    italy.NAME_1 = italy.NAME_1.replace({
        'Apulia': 'Puglia',
        'Emilia-Romagna': 'Emilia Romagna',
        'Friuli-Venezia Giulia': 'Friuli Venezia Giulia',
        'Sicily': 'Sicilia',
        'Trentino-Alto Adige': 'Trentino Alto Adige',
        'Valle d\'Aosta': 'Valle Aosta'
    })
    
    # Rename regions inside DFs and change column names
    for region, reg_df in regional_data.items():
        regional_data[region] = reg_df.rename(columns={
            'data': 'day',
            'denominazione_regione': 'region',
            'nuovi_tamponi': 'tested_swabs',
            'nuovi_positivi': 'positive_cases'
        })
        regional_data[region]['region'] = region
        
    return regional_data, pd.concat(regional_data.values())
#
###########################################################################

###########################################################################
# Utilties to create the instances for the model
def get_reg_day_demand(r, multiplier=1.00, overwhelmed_multiplier=1.20):
    """ Regional day demand.

        Corresponds to the real-life demand, up to a multiplier.
        For overwhelmed regions (determined based on news and press releases)
        there is a further, compound, multiplier.
    """

    sr = sorted(r.region.unique())
    ndays = len(r.day.unique())
    d = [list() for _ in sr]
    
    overwhelmed_regions = [
        8,  # Lombardia
        9,  # Marche
        11, # Piemonte
        19  # Veneto
    ]
    
    for idx, region in enumerate(sr):
        for day in range(ndays):
            abs_day = r.day.min() + timedelta(days=day)
            row = r[(r.region == region) & (r.day == abs_day)]
            
            if len(row) != 1:
                print(f"Error for region={region}, day={day}")
                
            demand = row.iloc[0].tested_swabs * multiplier
            
            if idx in overwhelmed_regions:
                demand *= overwhelmed_multiplier
                
            d[idx].append(int(demand))
    
    return d
            
def get_lab_day_demand(l, rdd):
    """ Daily demand, lab by lab.

        Assumes each region splits its demand equally to its labs.
    """

    sr = sorted(l.region.unique())
    ndays = len(rdd[0])
    d = [list() for _ in range(len(l))]
    
    for idx, labrow in l.iterrows():
        for day in range(ndays):
            reg_demand = rdd[sr.index(labrow.region)][day]
            reg_labs = len(l[l.region == labrow.region])
            d[idx].append(int(reg_demand / reg_labs))
            
    return d

def get_lab_closest_factory(l, f):
    """ Gives the factory closest to each lab. We assume labs procure reagents from
        their closest factory, at a rate mandated by the planner.
    """

    d = list()
    
    for idx, labrow in l.iterrows():
        lab_lon = [labrow.city_long] * len(f)
        lab_lat = [labrow.city_lat] * len(f)
        
        dist = [haversine(x, y, z, w) for x, y, z, w in zip(lab_lon, lab_lat, f.city_long.values, f.city_lat.values)]
        closest = dist.index(min(dist))
        
        d.append(closest)
        
    return d

def get_fac_served_labs(lcf):
    """ Inverse of get_lab_closest_factory. For each factory,
        gives the list of labs which have it as their closest.
    """

    nfac = max(lcf) + 1
    d = [list() for _ in range(nfac)]
    
    for facidx in range(nfac):
        for labidx, labfac in enumerate(lcf):
            if labfac == facidx:
                d[facidx].append(labidx)
                
    return d

def get_fac_day_demand(fsl, lr, l, r):
    """ Reagent demand for each factory, day by day.
        It sums the daily demands of all labs which procure from the factory.
    """

    nfac = len(fsl)
    ndays = len(r.day.unique())
    sr = sorted(r.region.unique())
    d = list()
        
    for faclabs in fsl:
        demand = [0] * ndays
        
        for day in range(ndays):
            abs_day = r.day.min() + timedelta(days=day)
            
            for lab in faclabs:
                region = sr[lr[lab]]
                rdem = r[(r.region == region) & (r.day == abs_day)].iloc[0].tested_swabs
                ldem = rdem / len(l[l.region == region])
                demand[day] += int(ldem)
        
        d.append(demand)
    
    return d

def get_fac_day_production(fdd, multiplier=1.2):
    """ Factory production, day by day. It's given by the demand,
        times a multiplier.
    """

    return [[int(x * multiplier) for x in row] for row in fdd]

def get_lab_region(l):
    """ List of labs in each region.
    """

    sr = sorted(l.region.unique())
    d = list()
    
    for _, labrow in l.iterrows():
        d.append(sr.index(labrow.region))
        
    return d

def get_lab_capacity(ldd, multiplier=1.25):
    """ Daily lab capacity. It is the maximum demand served by a lab,
        times a multiplier. If we know from real data that a lab has
        tested 5000 swabs in a day, then its capacity must be at least
        5000, so multiplier should be >= 1.0.
    """

    return [int(multiplier * max(labdemands)) for labdemands in ldd]

def get_lab_start_reagents(lc, multiplier=0.1):
    """ Initial number of reagents stored at each labs at the beginning
        of the planning horizon. If the planning horizon starts mid-emergency
        just set multiplier to 0, to have no safety stock.
    """
    return [int(multiplier * capacity) for capacity in lc]

def get_fac_start_reagents(fdd, multiplier=0.1):
    """ Initial number of reagents stored at each factory at the beginning
        of the planning horizon. If the planning horizon starts mid-emergency
        just set multiplier to 0, to have no safety stock.
    """

    return [int(multiplier * max(facdemands)) for facdemands in fdd]

def get_reg_max_inbound_reagents(rdd, multiplier=1.0):
    """ Number of reagents that can be delivered to each region every day,
        due to logistic constraints.
    """

    if multiplier == 'reg':
        return get_reg_max_inbound_reagents_realistic()
    else:
        return [int(multiplier * max(regdemands)) for regdemands in rdd]

def get_reg_max_inbound_swabs(rdd, multiplier=1.0):
    """ Number of swabs that can be delivered to each region every day,
        due to logistic constraints.
    """

    if multiplier == 'reg':
        return get_reg_max_inbound_swabs_realistic()
    else:
        return [int(multiplier * max(regdemands)) for regdemands in rdd]
    
def get_reg_max_inbound_reagents_realistic():
    """ Number of reagents that can be delivered to each region every day,
        due to logistic constraints. These figures are estimated taking into
        account geographical and economics characteristics of each region.
    """

    with_fac   = 10000 # If region has factory
    large_reg  =  4000 # If large/well-connected region
    normal_reg =  1000 # If normal region
    small_reg  =   500 # If small/impervious region

    return [
        small_reg,  # Abruzzo
        small_reg,  # Basilicata
        normal_reg, # Calabria
        large_reg,  # Campania
        with_fac,   # Emilia Romagna
        normal_reg, # Friuli Venezia Giulia
        with_fac,   # Lazio
        normal_reg, # Liguria
        with_fac,   # Lombardia
        normal_reg, # Marche
        small_reg,  # Molise
        with_fac,   # Piemonte
        large_reg,  # Puglia
        normal_reg, # Sardegna
        large_reg,  # Sicilia
        large_reg,  # Toscana
        small_reg,  # Trentino Alto Adige
        small_reg,  # Umbria
        small_reg,  # Valle Aosta
        with_fac    # Veneto
    ]

def get_reg_max_inbound_swabs_realistic():
    """ Number of swabs that can be delivered to each region every day,
        due to logistic constraints. These figures are estimated taking into
        account geographical and economics characteristics of each region.
    """

    large_reg  = 2000 # If large/well-connected region
    normal_reg = 100  # If normal region
    small_reg  = 500  # If small/impervious region
    
    return [
        small_reg,  # Abruzzo
        small_reg,  # Basilicata
        normal_reg, # Calabria
        large_reg,  # Campania
        large_reg,  # Emilia Romagna
        normal_reg, # Friuli Venezia Giulia
        large_reg,  # Lazio
        normal_reg, # Liguria
        large_reg,  # Lombardia
        normal_reg, # Marche
        small_reg,  # Molise
        large_reg,  # Piemonte
        large_reg,  # Puglia
        normal_reg, # Sardegna
        large_reg,  # Sicilia
        large_reg,  # Toscana
        small_reg,  # Trentino Alto Adige
        small_reg,  # Umbria
        small_reg,  # Valle Aosta
        large_reg   # Veneto
    ]

def get_lab_lab_compatible(l, inst_type):
    """ Lab-to-lab swab transfers compatibility map.
    """

    d = list()
    
    for idx1, labrow1 in l.iterrows():
        row = list()
        
        for idx2, labrow2 in l.iterrows():
            if idx1 == idx2:
                row.append(0)
                continue
                
            if labrow1.region == labrow2.region:
                row.append(1)
                continue
                
            if labrow1.region == 'Sardegna' or labrow2.region == 'Sardegna':
                row.append(0)
                continue
                
            if inst_type == 'regional':
                row.append(0)
            else:
                dist = haversine(labrow1.city_long, labrow1.city_lat, labrow2.city_long, labrow2.city_lat)
                
                if dist <= inst_type:
                    row.append(1)
                else:
                    row.append(0)
                
        d.append(row)
        
    return d

def get_fac_lab_compatible(nfac, lcf, inst_type):
    """ Factory-to-lab reagent delivery map.
    """

    d = list()
    
    for fac in range(nfac):
        if inst_type == 'regional':
            row = list()

            for lab, cfac in enumerate(lcf):
                if cfac == fac:
                    row.append(1)
                else:
                    row.append(0)
            
            d.append(row)
        else:
            d.append([1] * len(lcf))
            
    return d

def get_instance(f, l, r, inst_type, **kwargs):
    """ Generates an instance. It returns a dictionaty which, if
        dumped into a JSON file, can be read from the model solver.
    """

    lab_region = get_lab_region(l)
    
    real_data_reg_day_demand = get_reg_day_demand(r, multiplier=1.0)
    real_data_lab_day_demand = get_lab_day_demand(l, real_data_reg_day_demand)
    real_data_lab_capacity = get_lab_capacity(real_data_lab_day_demand, multiplier=1.0)
    
    reg_day_demand = get_reg_day_demand(r,
                                        multiplier=kwargs.get('reg_day_demand_mult', 1.0),
                                        overwhelmed_multiplier=kwargs.get('reg_overwhelmed_demand_mult', 1.0))
    lab_day_demand = get_lab_day_demand(l, reg_day_demand)
    
    lab_closest_factory = get_lab_closest_factory(l, f)
    fac_served_labs = get_fac_served_labs(lab_closest_factory)
    
    real_data_fac_day_demand = get_fac_day_demand(fac_served_labs, lab_region, l, r)
    lab_capacity = get_lab_capacity(real_data_lab_day_demand, multiplier=kwargs.get('lab_capacity_mult', 1.25))
    
    return dict(
        n_labs=len(l),
        n_factories=len(f),
        n_days=len(r.day.unique()),
        n_regions=len(r.region.unique()),
        lab_region=lab_region,
        lab_capacity=lab_capacity,
        lab_start_reagents=get_lab_start_reagents(real_data_lab_capacity, multiplier=kwargs.get('lab_start_reagents_mult', 0.1)),
        fac_start_reagents=get_fac_start_reagents(real_data_fac_day_demand, multiplier=kwargs.get('fac_start_reagents_mult', 0.1)),
        reg_max_inbound_reagents=get_reg_max_inbound_reagents(real_data_reg_day_demand, multiplier=kwargs.get('reg_max_in_reagents_mult', 1.0)),
        reg_max_inbound_swabs=get_reg_max_inbound_swabs(real_data_reg_day_demand, multiplier=kwargs.get('reg_max_in_swabs_mult', 1.0)),
        fac_day_production=get_fac_day_production(real_data_fac_day_demand, multiplier=kwargs.get('fac_day_production_mult', 1.2)),
        reg_day_demand=reg_day_demand,
        lab_lab_compatible=get_lab_lab_compatible(l, inst_type),
        fac_lab_compatible=get_fac_lab_compatible(len(f), lab_closest_factory, inst_type),
        reg_names=sorted(l.region.unique())
    )

def create_instance(f, l, r, inst_type, **kwargs):
    """ Creates an instance for a particular model (inst_type).
    """

    i = get_instance(f, l, r, inst_type, **kwargs)

    with open(f"{inst_path}/italy-{inst_type}.json", 'w') as file:
        file.write(json.dumps(i, cls=NpEncoder, indent=4))
        
    return i

def create_all_instances(f, l, r, **kwargs):
    """ Creates instances for all models.
    """

    for inst_type in ['regional', 100, 200, 400]:
        create_instance(f, l, r, inst_type, **kwargs)
#
###########################################################################

###########################################################################
# How to get f, l, r to pass to create_all_instances?
#
# f = pd.read_csv('italy_factories_coords.csv')
# l = pd.read_csv('italy_labs_coords.csv')
# rd, r = get_regional_data()
#
# Parameters used:
# reg_day_demand_mult=1.05,       
# reg_overwhelmed_demand_mult=2.0,
# lab_capacity_mult=1.25,         
# lab_start_reagents_mult=0.0,    
# fac_start_reagents_mult=0.0,    
# reg_max_in_reagents_mult='reg', 
# reg_max_in_swabs_mult='reg',    
# fac_day_production_mult=1.5
#
###########################################################################