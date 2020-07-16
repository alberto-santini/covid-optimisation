from scipy.stats import truncnorm
import numpy as np
import random
import json

def sqdist(pt1, pt2):
    return (pt1[0] - pt2[0])**2 + (pt1[1] - pt2[1])**2

def get_rnd_coords(n, low=0, high=100):
    return np.random.uniform(low=low, high=high, size=(n, 2)).tolist()

def get_rnd_coords_nearby(pt, radius=20, low=0, high=100):
    while True:
        theta = np.random.uniform(low=0, high=2*np.pi)
        r = np.sqrt(np.random.uniform(low=0, high=radius**2))
        x, y = pt[0] + r * np.cos(theta), pt[1] + r * np.sin(theta)
        
        if (low <= x <= high) and (low <= y <= high):
            return x, y
        
def get_labs(n_labs, regions, low=0, high=100):
    labs = list()
    
    for _ in range(n_labs):
        region = random.choice(regions)
        labs.append(get_rnd_coords_nearby(region, low=low, high=high))
        
    return labs

def get_lab_lab_distance(labs):
    dist = list()

    for l1 in labs:
        dist.append(list())

        for l2 in labs:
            d = np.sqrt(sqdist(l1, l2))
            dist[-1].append(np.round(d, decimals=2))

    return dist

def get_fac_lab_distance(labs, factories):
    dist = list()

    for f in factories:
        dist.append(list())

        for l in labs:
            d = np.sqrt(sqdist(f, l))
            dist[-1].append(np.round(d, decimals=2))

    return dist

def get_lab_regions(labs, regions):
    lab_regions = list()
    
    for lab in labs:
        closest = min(regions, key=lambda region: sqdist(lab, region))
        lab_regions.append(regions.index(closest))
        
    return lab_regions

def get_reg_demand(nr, nd, labs_per_region, base_demand, frac_peak_regions, peak_regions_mult, nonpeak_regions_mult, demand_noise_lower, demand_noise_upper):
    demand = [[base_demand] * nd for _ in range(nr)]
    
    for idx, drow in enumerate(demand):
        drow = [d * labs_per_region[idx] for d in drow]
        
        if random.random() <= frac_peak_regions:
            drow = [d * peak_regions_mult for d in drow]
        else:
            drow = [d * nonpeak_regions_mult for d in drow]
            
        drow = [d * truncnorm.rvs(demand_noise_lower, demand_noise_upper) for d in drow]
        drow = [int(d) for d in drow]
        
        demand[idx] = drow
        
    return demand

def get_lab_demand(nl, nd, reg_demand, lab_regions, labs_per_region):
    demand = [[0] * nd for _ in range(nl)]
    
    for idxl in range(nl):
        idxr = lab_regions[idxl]
                
        for idxd, dem in enumerate(reg_demand[idxr]):
            demand[idxl][idxd] = int(dem / labs_per_region[idxr])
            
    return demand

def get_labs_per_region(lr):
    lr = np.array(lr)
    nr = lr.max() + 1
    labs_per_region = []
    
    for region in range(nr):
        labs_per_region.append(np.array(lr[lr == region]).size)
        
    return labs_per_region

def get_lab_capacity(n, base_demand, lab_capacity_mult):
    lc = [int(base_demand * lab_capacity_mult) for _ in range(n)]
    return [x * truncnorm.rvs(0.75, 1.25) for x in lc]

def get_lab_fac_comp_l(labs, factories, lab_closest_facs):
    comp = list()
    
    for lab in labs:
        dists = [sqdist(lab, factory) for factory in factories]
        ordered_facs = np.array(dists).argsort()[:lab_closest_facs].tolist()
        comp.append(ordered_facs)
        
    return comp

def get_fac_lab_comp_l(n_factories, lfc):
    comp = list()
    
    for idxf in range(n_factories):
        clabs = [idxl for idxl, lcomp in enumerate(lfc) if idxf in lcomp]
        comp.append(clabs)
        
    return comp

def get_fac_lab_comp(n_labs, flc):
    return [[1 if idxl in row else 0 for idxl in range(n_labs)] for row in flc]

def get_lab_lab_comp(labs, lab_regions, lab_lab_comp_radius):
    nl = len(labs)
    comp = [[0] * nl for _ in labs]
    
    for idx1, lab1 in enumerate(labs):
        for idx2, lab2 in enumerate(labs):
            if idx1 == idx2:
                continue
                
            if lab_regions[idx1] == lab_regions[idx2]:
                comp[idx1][idx2] = 1
                continue
                
            if sqdist(lab1, lab2) <= lab_lab_comp_radius ** 2:
                comp[idx1][idx2] = 1
                continue
    
    return comp

def get_reg_capacity(n_regions, labs_per_region, reg_capacity_mult_swab_normal, reg_capacity_mult_reagent_normal, frac_hard_logistic_reg, reg_capacity_mult_swab_hard, reg_capacity_mult_reagent_hard, base_demand):
    sc = list()
    rc = list()
    
    for idxr in range(n_regions):
        smult = reg_capacity_mult_swab_normal
        rmult = reg_capacity_mult_reagent_normal
        
        if random.random() < frac_hard_logistic_reg:
            smult = reg_capacity_mult_swab_hard
            rmult = reg_capacity_mult_reagent_hard
            
        s = base_demand * labs_per_region[idxr] * smult
        r = base_demand * labs_per_region[idxr] * rmult
        
        sc.append(int(s))
        rc.append(int(r))
        
    return sc, rc

def get_lab_start_reagents(n, base_demand, lab_start_reagents_lower, lab_start_reagents_upper):
    s = [base_demand * truncnorm.rvs(lab_start_reagents_lower, lab_start_reagents_upper) for _ in range(n)]
    return[int(x) for x in s]

def get_fac_start_reagents(nf, nl, base_demand, fac_start_reagents_lower, fac_start_reagents_upper):
    s = [base_demand * nl / nf * truncnorm.rvs(fac_start_reagents_lower, fac_start_reagents_upper)]
    return [int(x) for x in s]

def get_fac_day_production(nf, nl, nd, base_demand, fac_prod_multiplier, fac_prod_lower, fac_prod_upper, fac_prod_type, n_release_days):
    p = list()
    
    for _ in range(nf):
        fp = list()
        
        for _ in range(nd):
            d = base_demand * nl / nf
            d *= fac_prod_multiplier
            d *= truncnorm.rvs(fac_prod_lower, fac_prod_upper)
            fp.append(int(d))
            
        if fac_prod_type == 'bumpy':
            release_days = random.sample(range(7), n_release_days)
            cum_demand = 0
            
            for day in range(nd):
                cum_demand += fp[day]
                dow = day % 7
                
                if dow in release_days:
                    fp[day] = cum_demand
                    cum_demand = 0
                else:
                    fp[day] = 0
                    
        p.append(fp)
        
    return p

if __name__ == "__main__":
    default_options = {
        'square_sz': 100,
        'square_boundary': 0.1,
        'n_labs': 100,
        'n_labs_per_region': 10,
        'n_facs_per_region': 0.5,
        'n_days': 14,
        'base_demand': 100,
        'frac_peak_regions': 0.40,
        'peak_regions_mult': 2.0,
        'nonpeak_regions_mult': 0.75,
        'demand_noise_lower': 0.95,
        'demand_noise_upper': 1.05,
        'lab_capacity_mult': 1.0,
        'lab_closest_facs': 1,
        'lab_lab_comp_radius': 10,
        'frac_hard_logistic_reg': 0.2,
        'reg_capacity_mult_reagent_normal': 1.0,
        'reg_capacity_mult_reagent_hard': 0.75,
        'reg_capacity_mult_swab_normal': 0.25,
        'reg_capacity_mult_swab_hard': 0.1,
        'fac_start_reagents_lower': 0.0,
        'fac_start_reagents_upper': 0.25,
        'lab_start_reagents_lower': 0.0,
        'lab_start_reagents_upper': 0.25,
        'fac_prod_multiplier': 1.0,
        'fac_prod_type': 'bumpy',
        'fac_bumpy_n_release_days': 2,
        'fac_prod_lower': 0.95,
        'fac_prod_upper': 1.05
    }

    inst_n = 1
    rows = list()
    low = 0
    high = 100
    options = default_options.copy()

    with open('inst-directory.csv', 'w') as d:
        strict_low = low + options['square_boundary'] * (high - low)
        strict_high = high - options['square_boundary'] * (high - low)
        n_labs = options['n_labs']
        
        for n_days in [5, 7, 10, 14]:
            for n_labs_per_region in [5, 10, 20]:
                n_regions = int(n_labs / n_labs_per_region)
                
                for n_facs_per_region in [0.1, 0.25, 0.5, 1]:
                    n_factories = max(1, int(n_regions * n_facs_per_region))
                    regions = get_rnd_coords(n_regions, low=strict_low, high=strict_high)
                    labs = get_labs(n_labs, regions, low=low, high=high)
                    factories = get_rnd_coords(n_factories, low=strict_low, high=strict_high)
                    lab_regions = get_lab_regions(labs, regions)
                    labs_per_region = get_labs_per_region(lab_regions)
                    lab_lab_distance = get_lab_lab_distance(labs)
                    fac_lab_distance = get_fac_lab_distance(labs, factories)
                    
                    reg_demand = get_reg_demand(
                        n_regions,
                        n_days,
                        labs_per_region,
                        options['base_demand'],
                        options['frac_peak_regions'],
                        options['peak_regions_mult'],
                        options['nonpeak_regions_mult'],
                        options['demand_noise_lower'],
                        options['demand_noise_upper'])
                    
                    lab_demand = get_lab_demand(n_labs, n_days, reg_demand, lab_regions, labs_per_region)

                    for lab_capacity_mult in [0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5]:
                        lab_capacity = get_lab_capacity(n_labs, options['base_demand'], lab_capacity_mult)
                        
                        for lab_closest_facs in [1, 2, 3]:
                            lab_fac_comp_l = get_lab_fac_comp_l(labs, factories, lab_closest_facs)
                            fac_lab_comp_l = get_fac_lab_comp_l(n_factories, lab_fac_comp_l)
                            fac_lab_comp = get_fac_lab_comp(n_labs, fac_lab_comp_l)
                            
                            for lab_lab_comp_radius in [0, 5, 10, 15, 20, 25]:
                                lab_lab_comp = get_lab_lab_comp(labs, lab_regions, lab_lab_comp_radius)
                                reg_capacity_swabs, reg_capacity_reagents = get_reg_capacity(
                                    n_regions,
                                    labs_per_region,
                                    options['reg_capacity_mult_swab_normal'],
                                    options['reg_capacity_mult_reagent_normal'],
                                    options['frac_hard_logistic_reg'],
                                    options['reg_capacity_mult_swab_hard'],
                                    options['reg_capacity_mult_reagent_hard'], 
                                    options['base_demand'])
                                lab_start_reagents = get_lab_start_reagents(
                                    n_labs, 
                                    options['base_demand'],
                                    options['lab_start_reagents_lower'],
                                    options['lab_start_reagents_upper'])
                                fac_start_reagents = get_fac_start_reagents(
                                    n_factories,
                                    n_labs,
                                    options['base_demand'],
                                    options['fac_start_reagents_lower'],
                                    options['fac_start_reagents_upper'])
                                
                                for fac_prod_multiplier in [0.8, 0.9, 1, 1.1, 1.2]:
                                    for fac_prod_type in ['steady', 'bumpy']:
                                        fac_day_production = get_fac_day_production(
                                            n_factories,
                                            n_labs,
                                            n_days,
                                            options['base_demand'],
                                            fac_prod_multiplier,
                                            options['fac_prod_lower'],
                                            options['fac_prod_upper'],
                                            fac_prod_type,
                                            options['fac_bumpy_n_release_days'])
                                        
                                        instance = {
                                            'n_labs': n_labs,
                                            'n_regions': n_regions,
                                            'n_factories': n_factories,
                                            'n_days': n_days,
                                            'lab_region': lab_regions,
                                            'lab_capacity': lab_capacity,
                                            'lab_start_reagents': lab_start_reagents,
                                            'fac_start_reagents': fac_start_reagents,
                                            'reg_max_inbound_reagents': reg_capacity_reagents,
                                            'reg_max_inbound_swabs': reg_capacity_swabs,
                                            'fac_day_production': fac_day_production,
                                            'reg_day_demand': reg_demand,
                                            'lab_lab_compatible': lab_lab_comp,
                                            'lab_lab_distance': lab_lab_distance,
                                            'fac_lab_compatible': fac_lab_comp,
                                            'fac_lab_distance': fac_lab_distance
                                        }

                                        with open(f"s-{inst_n}.json", 'w') as f:
                                            json.dump(instance, f)

                                        d.write(f"{inst_n},{n_days},{n_labs_per_region},{n_facs_per_region},")
                                        d.write(f"{lab_capacity_mult},{lab_closest_facs},{lab_lab_comp_radius},")
                                        d.write(f"{fac_prod_multiplier},{fac_prod_type}\n")
                                        
                                        inst_n += 1