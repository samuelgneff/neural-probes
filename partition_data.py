import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy.linalg as la
import scipy as sp


def segment_data(yws, yps, tolerance, deflection = 'large'):
#     print(yws.keys())
    elastic_limits_list = []
    slippage_limits_list = []
    ld_limits_list = []
    for key in yws.keys():
        print(yws.keys())
#         print('how many times')
        elastic_limits = RANSAC(yws['4'], yps['4'], 4, tolerance, deflection = 'large', graph=True)
        elastic_limits_list.append(elastic_limits)
#         print('elastic limits' ,elastic_limits)
        slippage_limits = find_slippage_segment(yws['4'], elastic_limits, threshold = 10)
        slippage_limits_list.append(slippage_limits)
#         print('slippage limits', slippage_limits)
        ld_limits = [elastic_limits[1], slippage_limits[0]]
        ld_limits_list.append(ld_limits)
#         print('ld limits', ld_limits )
        break
    return elastic_limits_list, slippage_limits_list, ld_limits_list



def RANSAC(xs, ys, post_num, tolerance = 2, deflection = 'large', graph = True):
    
    '''
    Calculates linear portion of small/large deflection curve
    
    '''
    if graph == True:
        
        fig, ax = plt.subplots(figsize=(15, 10))
    result = None
    distance = 0
    forward = True
    backward = True
    if deflection == 'small':
        
        start_point = len(ys)//2
        upper_limit = start_point + 2
        lower_limit = start_point - 2
        
    else:
        backward = False
        start_point = 0
        upper_limit = start_point + 2
        lower_limit = start_point
    
    y_sub = ys[lower_limit: upper_limit]
    x_sub = xs[lower_limit : upper_limit]
    for i in np.arange(0,len(ys)):
        
        
        if not backward and not forward:
            break
                
        if forward:
            upper_limit = start_point + i + 2
            if upper_limit > len(xs):
                forward = False
                continue
            
            # consider removing try catch statement if statement above may be 
            # sufficient.
            try:
                y_sub = ys[lower_limit: upper_limit]
                x_sub = xs[lower_limit: upper_limit]
            except IndexError:
                y_sub = ys[lower_limit: upper_limit - 1]
                x_sub = xs[lower_limit: upper_limit - 1] 
                forward = False
                continue
            
            line_obj = np.poly1d(np.polyfit(x_sub, y_sub, 1))
            result = line_obj(x_sub)
            p1 = np.array([x_sub[0], y_sub[0]])
            p2 = np.array([x_sub[-2], y_sub[-2]])
            p3 = np.array([x_sub[-1], y_sub[-1]])
            distance = (np.abs(np.cross(p2-p1, p3-p1)) / la.norm(p2-p1))
            adjusted_tolerance = tolerance  / (len(y_sub))
            if distance > adjusted_tolerance:
                forward = False
                distance = 0
                continue
            if graph == True:
                plotvr(fig, ax, result, x_sub, y_sub, xs, ys, post_num, i)
        if backward:
            lower_limit = start_point  - i - 1
            
            if lower_limit < 0:
                backward = False
                continue  
            try:
                y_sub = ys[lower_limit: upper_limit]
                x_sub = xs[lower_limit: upper_limit]
            except IndexError:
                y_sub = ys[lower_limit + 1: upper_limit]
                x_sub = xs[lower_limit + 1: upper_limit]
                backward = False
                continue
            line_obj = np.poly1d(np.polyfit(x_sub, y_sub, 1))
            result = line_obj(x_sub)
            p1 = np.array([x_sub[0],  y_sub[0]])
            p2 = np.array([x_sub[-2], y_sub[-2]])
            p3 = np.array([x_sub[1], y_sub[1]])
            distance = (np.abs(np.cross(p2-p1, p3-p1)) / la.norm(p2-p1))
            if distance > tolerance:
                backward = False
                distance = 0
                continue
            if graph == True:  
                plotvr(fig, ax, result, x_sub, y_sub, xs, ys, post_num)         
    return np.array([lower_limit, upper_limit])