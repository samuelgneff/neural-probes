import numpy.linalg as la
import numpy as np

def dual_deflection_analysis(folder = '186_Uninfiltrated_Row 7'):
    # grab all data files
    wire_files, ref_files, filenames = upload_tracker_files(folder)
    # find all yws and yps
    yps, yws, post_numbers = find_yps_yws(filenames, wire_files, ref_files)
    
    # segment data 
    elastic_limits, slippage_limits, large_deflection_limits = segment_data(yws, yps, tolerance = 8, deflection = 'large')
    
    #plot_basic(yps,yws)
    #plot_final_segmentation(yws, yps, elastic_limits, slippage_limits, large_deflection_limits, )

    #LD_model(yws['12'], yps['12'])
    
    # calculate TLS fit of large deflection
    #y_calc, x_error, a_tls, fro_norm = tls(total_distance_dict['wire15'], yp)
    
    
    
    
    return elastic_limits , slippage_limits, large_deflection_limits, yps, yws


def tls(X,y):
    
    if X.ndim == 1: 
        n = 1 # the number of variable of X
        X = X.reshape(len(X),1)
    else:
        n = np.array(X).shape[1] 
    
    Z = np.vstack((X.T,y)).T
    U, s, Vt = la.svd(Z, full_matrices=True)
    
    V = Vt.T
    Vxy = V[:n,n:]
    Vyy = V[n:,n:]
    a_tls = - Vxy  / Vyy # total least squares soln
    
    Xtyt = - Z.dot(V[:,n:]).dot(V[:,n:].T)
    Xt = Xtyt[:,:n] # X error
    y_tls = (X+Xt).dot(a_tls)
    fro_norm = la.norm(Xtyt, 'fro')
    
    return y_tls, X+Xt, a_tls, fro_norm
    

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



def LD_model(YW,YP):
    # ===================== calculations =======================

    #Constants
    gamma = 0.85
    Kt = 2.65

    #wire and post contstants
    Ew = 213e9#Pa
    Lw = 3900 # length of wire in um
    rw = 19/2 # radius of wire in um
    Iw = (np.pi/4)*rw**4

    Ep = 0.621e6
    Lp = 270 # um - difference of total length (536 um) and offset (30 um)
    rp = 39/2
    Ip = (np.pi/4)*rp**4
    yp = np.linspace(0,220,100)

    n = 0
    phi = np.pi / 2

    #some parameter definitions
    eta = np.sqrt(1 + n**2)
    theta = np.arcsin(yp / (gamma*Lp))

    c = Kt / (3*eta)
    A = Ip / (Ew*Iw)
    B = (Lw**3) / (Lp**2)
    C = theta / (np.sin(phi - theta))

    yw = c * Ep * A * B * C


    #================ Calcualte Ep (from post 15) ==================

#     YW = YWYP(:,1);
#     YP = YWYP(:,2);

    eta = np.sqrt(1 + n**2)
    theta = np.arcsin(YP / (gamma*Lp))

    c = Kt / (3*eta)
    A = Ip / (Ew*Iw)
    B = (Lw**3) / (Lp**2)
    C = theta / (np.sin(phi - theta))

    EP = YW / (c * A * B * C) #Post modulud calculated at each deflection position

    #================= Uncertainty Quanitification ================ 

    N = 1000 # sample size for monte carlo

    Ew_sd = Ew*0.005
    Lw_sd = 100/2
    rw_sd = 1.27/2

    Ep_sd = Ep*0.005
    Lp_sd = 5/2
    rp_sd = 1/2
    # yp_sd = linspace(0,)


    # Monte Carlo Sampling
    Ew_mc = np.random.normal(Ew,Ew_sd,(N,1)) # random normally distributed sampling of values
    Lw_mc = np.random.normal(Lw,Lw_sd,(N,1))
    rw_mc = np.random.normal(rw,rw_sd,(N,1))
    Iw_mc = (np.pi/4)*np.power(rw_mc,4) # Moment Inertia of circular cross-section

    Ep_mc = np.random.normal(Ep,Ep_sd,(N,1)) # Same as above for the post
    Lp_mc = np.random.normal(Lp,Lp_sd,(N,1))
    rp_mc = np.random.normal(rp,rp_sd,(N,1))
    Ip_mc = (np.pi/4)*np.power(rp_mc,4)

    yp_mc = np.linspace(0,210,1000); # Evenly spaced deflection from 0 to 210 microns (1000 points)

    n = 0  #no lateral force, just vertical
    phi = np.pi / 2  # initial angle

    # some parameter definitions
    eta = np.sqrt(1 + n**2)
    theta = np.arcsin(yp_mc / (gamma*Lp_mc))

    # Derived equation for large dflection relationship between yw and yp
    c = Kt / (3*eta)
    A = Ip_mc / (Ew_mc*Iw_mc)
    B = (np.power(Lw_mc,3)) / np.power(Lp_mc,2)
    C = theta / (np.sin(phi - theta))

    yw_mc = c * Ep_mc * A * B * C
    [m,n] = yw_mc.shape
    fig, ax = plt.subplots(figsize=(10, 10))

    for i in  np.arange(0,len(yp_mc)):
        ax.plot(np.ones((m,1))*yp_mc[i],yw_mc[:,i], color='blue', zorder=1)
    y_calc, x_error, a_tls, fro_norm = tls(yw, yp)
    
    ax.plot(y_calc, yw, label='TLS fit to Modulus', color='green', zorder=2)
    ax.scatter(yp, yw, label='data', color='red', zorder=2)
    ax.set_xlabel('x axis', fontsize = '20')
    ax.set_ylabel('y axis', fontsize = '20')
    ax.legend()
    plt.show()