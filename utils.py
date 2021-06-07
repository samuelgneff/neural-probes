def remove_offset(initial_values, data):

    ''' A function to remove the offset from Tracker
    
    '''


    return_data = data.copy()
    for i in range(len(data)):
        return_data[i][1] = (initial_values[0] - data[i][1])
        return_data[i][2] = (initial_values[1] - data[i][2])
    return return_data


def plot_final_segmentation(yws, yps, elastic_lim, slippage_lim, ld_lim): #, max_val, max_x, min_val, min_x, Fc):
    
    '''
    Plot utility function. 
    
    '''
    for i, key in enumerate(yws.keys()):
        yp = yps[key]
        yw = yws[key]
        fig, ax = plt.subplots(figsize=(15, 10))
        ax.set_title('Segmented Deflection Curve, post - {0}'.format(key), fontsize = '20')
        ax.scatter(yp[elastic_lim[i][0]:elastic_lim[i][1]], yw[elastic_lim[i][0]:elastic_lim[i][1]], label = 'elastic', color='green')
        ax.scatter(yp[ld_lim[i][0]:ld_lim[i][1]], yw[ld_lim[i][0]:ld_lim[i][1]], label = 'large deflection', color='orange')
        ax.scatter(yp[slippage_lim[i][0]:slippage_lim[i][1]], yw[slippage_lim[i][0]:slippage_lim[i][1]], label = 'slippage')
    #     ax.plot(result,x_sub , c = 'green', lw = '3')
    #     ax.plot(y_sub[0],x_sub[0] , 'or')
    #     ax.plot(y_sub[-1], x_sub[-1], 'oy')
        ax.set_xlabel('yp (microns)', fontsize = '20')
        ax.set_ylabel('yw (microns)', fontsize = '20')
        ax.legend()
        fig.savefig('Segmented Deflection Curve, post - {0}'.format(key))
        ax.plot()


def plot_basic(yws, yps): #, max_val, max_x, min_val, min_x, Fc):
    
    '''
    Plot utility function. 
    
    '''
    for i, key in enumerate(yws.keys()):
        yp = yps[key]
        yw = yws[key]
        fig, ax = plt.subplots(figsize=(15, 10))
        ax.set_title('Deflection Curve, post - {0}'.format(key), fontsize = '20')
        ax.scatter(yw,yp, label='Deflection')
        ax.set_xlabel('yp (microns)', fontsize = '20')
        ax.set_ylabel('yw (microns)', fontsize = '20')
        ax.legend()
        fig.savefig('Deflection Curve, post - {0}'.format(key))
        ax.plot()


def find_slippage_segment(wireData, elastic_limits, threshold = 10, graph=True):
    slippage_start = elastic_limits[1]
    slippage_end = len(wireData)
#     print('len(wireData) ', len(wireData))
#     print('wire data', wireData)
    
    for i in range(len(wireData[elastic_limits[1]:])):
        if i == len(wireData[elastic_limits[1]:]) -1:
            break
        wireStepSize = abs(wireData[elastic_limits[1] + (i+1)] - wireData[elastic_limits[1] + (i)])
        
#         print('wirestepsize: ', wireStepSize)
#         print('indice: ', elastic_limits[1] + (i+1))
        #plotslippagecalc(fig, ax, x_sub, y_sub, xs, ys, post_num)
        if (wireStepSize > threshold):
#             print('do we enter this if?')
            slippage_start = elastic_limits[1] + i + 1
            slippage_end = len(wireData)
#             wireData = wireData[:-i]
            break
        #wireData, refData = makeSameLength(wireData, refData)
    return [slippage_start, slippage_end]

def find_yp_yw(wire, ref):
    #     print(len(wire))
    #     print(len(ref))
    yw = process_data(wire)
    ref = process_data(ref)
    yp = abs(np.subtract(yw,  ref))

    return yw, yp


def upload_tracker_files(folder):
    filenames = sorted(os.listdir(folder))
    wire_files_data = []
    ref_files_data = []
    for i in np.arange(0, len(filenames), 2):
        #print(filenames[i])
        #print(filenames[i+1])
        wire_data = np.loadtxt('{0}/{1}'.format(folder, filenames[i+1]), skiprows = 2)
        wire_data_ref = np.loadtxt('{0}/{1}'.format(folder, filenames[i]), skiprows =2)
        wire_data, wire_data_ref = make_same_length(wire_data, wire_data_ref)
#         print('len of wire_data: ', len(wire_data))
#         print('len of wire_data_ref: ', len(wire_data_ref))
#         print(filenames[i])
#         print(filenames[i+1])
        wire_files_data.append(wire_data)
        ref_files_data.append(wire_data_ref) 
        
    print('length of wire files: ', len(wire_files_data))
    print('length of ref files: ', len(ref_files_data))
    return wire_files_data, ref_files_data, filenames