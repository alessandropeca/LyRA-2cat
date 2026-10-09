
def LYRA_LR_compl_part4(args):
    globals().update(args)
    global save_output
    global save_images
    global path_output
    global filename_outinfo
    global r_in_tm
    global r_min_nm
    global r_max_nm
    global r_lr
    global ra_input
    global dec_input
    global ra_output
    global dec_output
    global sigma_input
    global sigma_output
    global input_cat_type
    global ID_input
    global mags_input
    global lr_outID
    global lr_outRA
    global lr_outDEC
    global lr_inID
    global lr_inRA
    global lr_inDEC
    global lr_mag_in
    global lr_r
    global bfactor
    global prior
    global lr_LR
    global LR_threshold
    global rel
    global Nout
    global Nout_int
    global save_images
    global path_images
    global needed_only_once
    global e
    global f
    global lr_noIDerr
    global sigma_match
    global lr_noID_sec
    global origin
    global goto_pos_match
    global outID_names
    global inID_names

    import numpy as np
    import matplotlib.pyplot as plt
    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib
    matplotlib.use('Agg')
    # matplotlib.use('TkAgg')
    from scipy.interpolate import interp1d
    import IPython

    from LYRA_LR_Re_part3 import LYRA_LR_Re_part3
    globals().update(LYRA_LR_Re_part3(globals()))
    
    if goto_pos_match == False:
        ############################################################# Completness ###
        print('\nComputing thresholds:')
        if np.count_nonzero(lr_mag_in != nomag) >= min_sample_size:
            distr_quantile = lr_LR[np.logical_not(np.isnan(lr_LR))]
            quantile1 = np.quantile(distr_quantile, 0.01)
            quantile2 = np.quantile(distr_quantile, 0.99)
            LRth = np.linspace(quantile1, quantile2, 100)
            
            RC = []
            R = []
            C = []
            for i in range(len(LRth)):
                compl_rel = []
                compl_id = []
                for j in range(len(lr_LR)):
                    if lr_LR[j] > LRth[i]:
                        compl_rel.append(rel[j])
                        compl_id.append(lr_inID[j])
                if len(np.unique(compl_id)) == 0:
                    R_ = 0
                else:
                    R_ = np.sum(compl_rel)/len(np.unique(compl_id))
                R.append(R_)
                if Nout_int == 0:
                    C_ = 0
                else:
                    C_ = np.sum(compl_rel)/Nout_int
                C.append(C_)
                RC.append((R_+C_)/2)

            RC = np.where(np.isnan(RC), 0, RC)
            RCmax_index = np.where(RC==np.max(RC))
            RCmax_index = RCmax_index[0]
            
            if len(RCmax_index) > 1:
                RCmax_index = int(RCmax_index[0])
            else:
                RCmax_index = int(RCmax_index)
    
            if custom_RC_threshold != -99:
                max_value, max_index = lyra.find_plateau(LRth, RC, threshold=custom_RC_threshold)
                RCmax_index = max_index
            
            LR_threshold = LRth[RCmax_index]
            # # LR_threshold = LRth[2] 
            # file_prim = file_primary.split("/")[-1]
            # post_thresholds = {
            #                     'subaru_HSC_S82.csv': {
            #                         'g_kronflux_mag': 0.34,
            #                         'r_kronflux_mag': 0.32,
            #                         'i_kronflux_mag': 0.33,
            #                         'z_kronflux_mag': 0.33,
            #                         'y_kronflux_mag': 0.41,
            #                     },
            #                     'iacstr82_Fliri+16.csv': {
            #                         'MAG_AUTO_u': 0.05,
            #                         'MAG_AUTO_g': 0.04,
            #                         'MAG_AUTO_r': 0.04,
            #                         'MAG_AUTO_i': 0.04,
            #                         'MAG_AUTO_z': 0.04,
            #                     },
            #                     'SDSS_coadd_Jiang+14_combined_rgizu.csv': {
            #                         'MAG_AUTO_u_': 0.05,
            #                         'MAG_AUTO_g_': 0.03,
            #                         'MAG_AUTO_r_': 0.03,
            #                         'MAG_AUTO_i_': 0.03,
            #                         'MAG_AUTO_z_': 0.03,
            #                     },
            #                     'VHS_McMahon+20_2.csv': {
            #                         'Yapc4': 0.05,
            #                         'Japc4': 0.03,
            #                         'Hapc4': 0.02,
            #                         'Ksapc4': 0.03
            #                     },
            #                     'UKIDSS_DR11.csv': {
            #                         'YAPERMAG4': 0.06,
            #                         'JAPERMAG4': 0.03,
            #                         'HAPERMAG4': 0.03,
            #                         'KAPERMAG4': 0.03
            #                     },
            #                     'spiesall_timlin+16.csv': {
            #                         'mag_auto_ch1': 0.11,
            #                         'mag_auto_ch2': 0.16,
            #                     },
            #                     'shela_irac_v1.3_Papovich+16.csv': {
            #                         'mag3p6_auto': 0.11,
            #                         'mag4p5_auto': 0.15,
            #                     },
            #                     'AllWISE_S82_2.csv': {
            #                         'w1mproAB': 0.24,
            #                         'w2mproAB': 0.26,
            #                     },
            #                     'GALEX_GR67_FULL.csv': {
            #                         'FUV_MAG_AUTO': 0.05,
            #                         'NUV_MAG_AUTO': 0.05,
            #                     },
            #                     'herschell_viero+14.csv': {
            #                         'log_F250': 0.05,
            #                         'log_F350': 0.05,
            #                     },
            #                     'GAIA_DR3_S82_2.csv': {
            #                         'phot_g_mean_magAB': 0.05,
            #                         'phot_bp_mean_magAB': 0.05,
            #                         'phot_rp_mean_magAB': 0.05,
            #                     },
            #                     'FIRST_S82.csv': {
            #                         'log_FINT': 0.05
            #                     }
            # }
            # # LR_threshold = 0.05
            # LR_threshold = post_thresholds.get(file_prim, {}).get(smag, 0.05)
                
            # # intersection
            # diff = np.array(R) - np.array(C)
            # sign_changes = np.where(np.diff(np.sign(diff)))[0]
            # intersections = []
            # for index in sign_changes:
            #     x_vals = LRth[index:index+2]
            #     y_vals = diff[index:index+2]

            #     # Create a linear interpolation function
            #     f = interp1d(y_vals, x_vals)

            #     # Use the interpolation function to find the x corresponding to y=0
            #     intersect_x = f(0)
            #     intersections.append(intersect_x)

            # # IPython.embed()
            # LR_threshold = intersections[0]
            ###
            
            print('Pth: ' + str('{:.3g}'.format(LRth[RCmax_index])))
            print('C[Pth]: ' + str('{:.3g}'.format(C[RCmax_index])))
            print('R[Pth]: ' + str('{:.3g}'.format(R[RCmax_index])))

            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, '(R+C)/2 max: ' + str(round(RC[RCmax_index],2)) + '\n'
                            + 'Pth: ' + str('{:.3g}'.format(LRth[RCmax_index])) + '\n'
                            + 'C[Pth]: ' + str('{:.3g}'.format(C[RCmax_index])) + '\n'
                            + 'R[Pth]: ' + str('{:.3g}'.format(R[RCmax_index])))

            if save_images == True:
                
                selected_x = None
                selected_y = None
                # Function to capture and print the selected point
                def onpick(event):
                    global selected_x, selected_y
                    selected_x = event.mouseevent.xdata
                    selected_y = event.mouseevent.ydata

                fig6, ax = plt.subplots()
                # fig6.canvas.mpl_connect('pick_event', onpick)
                # plt.gca().set_picker(True)
                ax.plot(LRth, RC, color='r', label='(R+C)/2', zorder=3)
                ax.plot(LRth, R, color='mediumblue', ls='--', label='R', zorder=3)
                ax.plot(LRth, C, color='mediumblue', ls=':', label='C', zorder=3)
                ax.axvline(color='r', lw=1, ls='-.', alpha=0.8, x=LR_threshold)
                ax.set_xlabel('$LR_{th}$')
                ax.grid(ls=':', color='grey', alpha=0.3, zorder=0)
                ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
                ax.set_title('Completeness ' + smag)
                # ax.set_xscale('log')
                ax.legend()
                # plt.show()
                # Print the selected coordinates outside the function
                # if selected_x is not None and selected_y is not None:
                #     print(f'Selected point: X = {selected_x}, Y = {selected_y}')
                # else:
                #     print('No point selected.')
                fig6.savefig(path_images+'completeness_'+smag+add_name_output+'.png',bbox_inches="tight", dpi=200)

            # from indeces to names:
            lr_outID = lr_outID.astype(np.int32)
            lr_inID = lr_inID.astype(np.int32)
            outID_names = ID_output_names[lr_outID]
            inID_names = ID_input_names[lr_inID]
            outID_names = outID_names.astype(object)
            inID_names = inID_names.astype(object)
            # add non-matched sources
            lr_outID = np.append(lr_outID, lr_noID)
            lr_inID = np.append(lr_inID, np.full((len(lr_noID),),-99))
            outID_names = np.append(outID_names, lr_noID)
            inID_names = np.append(inID_names, np.full((len(lr_noID),),-99))
            
            lr_outRA = np.append(lr_outRA, lr_noIDra)
            lr_outDEC = np.append(lr_outDEC, lr_noIDdec)
            sigma_match = np.append(sigma_match, lr_noIDerr)
            # print("origin", origin.shape, type(origin))
            # print("lr_noID_sec", len(lr_noID_sec), type(lr_noID_sec))
            origin = list(origin)
            origin.extend(lr_noID_sec)
            # print("origin", origin.shape, type(origin))
            lr_inRA = np.append(lr_inRA, np.full((len(lr_noID),),-99))
            lr_inDEC = np.append(lr_inDEC, np.full((len(lr_noID),),-99))
            lr_mag_in = np.append(lr_mag_in, np.full((len(lr_noID),),-99))
            lr_r = np.append(lr_r, np.full((len(lr_noID),),-99))
            bfactor = np.append(bfactor, np.full((len(lr_noID),),-99))
            prior = np.append(prior, np.full((len(lr_noID),),-99))
            lr_LR = np.append(lr_LR, np.full((len(lr_noID),),-99))
            rel = np.append(rel, np.full((len(lr_noID),),-99))
        
        else:
            print("\n### Warning: Not enough good values for mag_input: " + smag + " inside the matching radius, so I skip it. ###")
            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, "\n### Warning: Not enough good values for mag_input: " + smag + " inside the matching radius, so I skip it. ###")
            
            # fill with -99s: e.g., bfactor is now only Bpos, which is already included in the calcs if Bphot is computed. 
            # If Bphot is never computed, the code goes to bayesyan pos. match, which computes again Bpos.
            LR_threshold = nomag
            bfactor[:] = -99
            prior[:] = -99
            lr_LR[:] = -99
            rel = np.array(rel)
            rel[:] = -99

            # from indeces to names:
            lr_outID = lr_outID.astype(np.int32)
            lr_inID = lr_inID.astype(np.int32)
            outID_names = ID_output_names[lr_outID]
            inID_names = ID_input_names[lr_inID]
            # add non-matched sources
            lr_outID = np.append(lr_outID, lr_noID)
            outID_names = np.append(outID_names, lr_noID)
            outID_names = outID_names.astype(object)
            lr_inID = np.append(lr_inID, np.full((len(lr_noID),),-99))
            inID_names = np.append(inID_names, np.full((len(lr_noID),),-99))
            inID_names = inID_names.astype(object)

            lr_outRA = np.append(lr_outRA, lr_noIDra)
            lr_outDEC = np.append(lr_outDEC, lr_noIDdec)
            sigma_match = np.append(sigma_match, lr_noIDerr)
            origin.extend(lr_noID_sec)
            lr_inRA = np.append(lr_inRA, np.full((len(lr_noID),),-99))
            lr_inDEC = np.append(lr_inDEC, np.full((len(lr_noID),),-99))
            lr_mag_in = np.append(lr_mag_in, np.full((len(lr_noID),),-99))
            lr_r = np.append(lr_r, np.full((len(lr_noID),),-99))
            bfactor = np.append(bfactor, np.full((len(lr_noID),),-99))
            prior = np.append(prior, np.full((len(lr_noID),),-99))
            lr_LR = np.append(lr_LR, np.full((len(lr_noID),),-99))
            rel = np.append(rel, np.full((len(lr_noID),),-99))

    return globals()

