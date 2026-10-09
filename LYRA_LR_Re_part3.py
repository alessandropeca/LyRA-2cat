
def LYRA_LR_Re_part3(args):
    globals().update(args)
    global r_out_tree
    global r_in_tree
    global r_lr
    global save_output
    global save_images
    global path_output
    global filename_outinfo
    global r_in_tm
    global r_min_nm
    global r_max_nm
    global r_lr
    global nm_interp
    global qm_interp
    global ra_input
    global dec_input
    global ra_output
    global dec_output
    global sigma_input
    global sigma_output
    global input_cat_type
    global ID_input
    global mags_input
    global Q_int
    global lr_outID
    global lr_outRA
    global lr_outDEC
    global lr_inID
    global lr_inRA
    global lr_inDEC
    global lr_mag_in
    global lr_r
    global lr_fr
    global qm_
    global nm_
    global bfactor
    global prior
    global lr_LR
    global rel
    global lr_noID
    global lr_noIDra
    global lr_noIDdec
    global needed_only_once
    global bayes
    global area_input
    global area_output
    global add_prior
    global qn_ratio_f
    global lr_noIDerr
    global sigma_match
    global origin
    global lr_noID_sec
    global d
    global goto_pos_match

    import numpy as np
    import matplotlib.pyplot as plt
    import warnings

    import LYRA_functions as lyra
    from LYRA_distrib_part2 import LYRA_distrib_part2
    globals().update(LYRA_distrib_part2(globals()))
    global tree
    import matplotlib
    matplotlib.use('Agg')
    
    warnings.filterwarnings("ignore")
    if goto_pos_match == False:
        ############################################################# Likelihood and Reliability ###
        print('\nLikelihood and Reliability:')
        neigh_idx_lr = tree.query_ball_point(r_out_tree, r=r_lr, workers=-1)
        len_neigh_idx_lr = tree.query_ball_point(r_out_tree, r=r_lr, workers=-1, return_length=True)
        lr_sum = np.sum(len_neigh_idx_lr)
        
        lr_inRA = np.zeros(lr_sum)
        lr_inDEC = np.zeros(lr_sum)
        lr_outRA = np.zeros(lr_sum)
        lr_outDEC = np.zeros(lr_sum)
        lr_mag_in = np.zeros(lr_sum)
        bfactor = np.zeros(lr_sum)
        prior = np.zeros(lr_sum)
        lr_LR = np.zeros(lr_sum)
        lr_inID = np.zeros(lr_sum)
        lr_outID = np.zeros(lr_sum)
        lr_r = np.zeros(lr_sum)
        lr_fr = np.zeros(lr_sum)
        qn_ratio = np.zeros(lr_sum)
        sigma_match = np.zeros(lr_sum)
        # origin = np.empty(lr_sum, dtype="<U9999")
        # origin = np.zeros(lr_sum, dtype="<S9999")
        origin = [0] * lr_sum
        rel = []
        lr_noID = []
        lr_noIDra = []
        lr_noIDdec = []
        lr_noIDerr = []
        lr_noID_sec = []

        d=0
        for i in range(len(neigh_idx_lr)):
            lr_single_outID = []
            for j in range(len(neigh_idx_lr[i])):
                rain_tmp = r_in_tree[neigh_idx_lr[i],0]
                decin_tmp = r_in_tree[neigh_idx_lr[i],1]
                Dra_tmp = (r_out_tree[i,0] - rain_tmp[j])
                Ddec_tmp = (r_out_tree[i,1] - decin_tmp[j])
                rr_lr = lyra.quadratic_sum(Dra_tmp, Ddec_tmp)
                if rr_lr <= r_lr: # just a check, the if condition should be already verified
                    ra_tmp = ra_input[neigh_idx_lr[i]]
                    lr_inRA[d] = ra_tmp[j]
                    dec_tmp = dec_input[neigh_idx_lr[i]]
                    lr_inDEC[d] = dec_tmp[j]
                    mag_tmp = mag_input[neigh_idx_lr[i]]
                    lr_mag_in[d] = mag_tmp[j]
                    inID_tmp = ID_input[neigh_idx_lr[i]]
                    lr_inID[d] = inID_tmp[j]
                    lr_outID[d] = ID_output[i]
                    lr_outRA[d] = ra_output[i]
                    lr_outDEC[d] = dec_output[i]
                    lr_r[d] = rr_lr
                    if input_cat_type == 1:
                        sigma_input_tmp = sigma_input[j]
                    else:
                        sigma_input_tmp = sigma_input
                    if output_cat_type == 1:
                        sigma_output_tmp = sigma_output[i]
                    else:
                        sigma_output_tmp = sigma_output
                    sigma_match[d] = min(sigma_output_tmp, sigma_input_tmp)
                    # pri_origin_tmp = pri_origin[neigh_idx_lr[i]]
                    pri_origin_tmp = [pri_origin[i_] for i_ in neigh_idx_lr[i]]
                    
                    if isinstance(pri_origin_tmp[j], tuple):
                        if isinstance(sec_origin[i], tuple):
                            origin[d] = pri_origin_tmp[j] + sec_origin[i] # merging two tuples
                        else:
                            if isinstance(sec_origin[i], str):
                                new_tup_sec = tuple(int(x) for x in sec_origin[i].strip("()").split(",") if x)
                                origin[d] = pri_origin_tmp[j] + new_tup_sec # merging two tuples
                            else:
                                print("\n\nAAA ERROR: sec_origin[i] is not a tuple:", sec_origin[i], type(sec_origin[i]))
                                sys.exit()
                    else:
                        if isinstance(pri_origin_tmp[j], str):
                            new_tup_pri = tuple(int(x) for x in pri_origin_tmp[j].strip("()").split(",") if x)
                            if isinstance(sec_origin[i], tuple):
                                origin[d] = new_tup_pri + sec_origin[i] # merging two tuples
                            elif isinstance(sec_origin[i], str):
                                new_tup_sec = tuple(int(x) for x in sec_origin[i].strip("()").split(",") if x)
                                origin[d] = new_tup_pri + new_tup_sec # merging two tuples
                            else:
                                print("\n\nAAA ERROR: sec_origin[i] is not a tuple:", sec_origin[i], type(sec_origin[i]))
                                sys.exit()		
                        else:
                            print("\n\nAAA ERROR: pri_origin_tmp[j] is not a tuple:", pri_origin_tmp[j], type(pri_origin_tmp[j]))
                            sys.exit()

                    # pri_origin_tmp_str = pri_origin_tmp[j]
                    # if isinstance(pri_origin_tmp_str, str) == False:
                    # 	pri_origin_tmp_str = str(pri_origin_tmp_str)
                    # sec_origin_str = sec_origin[i]
                    # if isinstance(sec_origin_str, str) == False:
                    # 	sec_origin_str = str(sec_origin_str)
                    # underscore_str = "_"
                    # pri_origin_tmp_str.replace('lyra_', '').replace('.csv', '')
                    # sec_origin_str.replace('lyra_', '').replace('.csv', '')
                    # origin[d] = (pri_origin_tmp_str + underscore_str + sec_origin_str).encode()
                    
                    if bayes == False: # standard LR calculation
                        if (input_cat_type == 2) or (input_cat_type == 0):
                            if (output_cat_type == 2) or (output_cat_type == 0):
                                lr_fr[d] = lyra.FR(sigma_output,sigma_input,Dra_tmp,Ddec_tmp)
                            else:
                                lr_fr[d] = lyra.FR(sigma_output[i],sigma_input,Dra_tmp,Ddec_tmp)
                            if lr_mag_in[d] != nomag:
                                qn_ratio[d] = qn_ratio_f(lr_mag_in[d])
                                lr_LR[d] = lyra.LR(lr_fr[d],qn_ratio[d])
                            else:
                                qn_ratio[d] = 0
                                lr_LR[d] = 0
                        
                        elif input_cat_type == 1:
                            if (output_cat_type == 2) or (output_cat_type == 0):
                                sigma_in_tmp = sigma_input[neigh_idx_lr[i]]
                                lr_fr[d] = lyra.FR(sigma_output,sigma_in_tmp[j],Dra_tmp,Ddec_tmp)
                            else:
                                sigma_in_tmp = sigma_input[neigh_idx_lr[i]]
                                lr_fr[d] = lyra.FR(sigma_output[i],sigma_in_tmp[j],Dra_tmp,Ddec_tmp)
                            if lr_mag_in[d] != nomag:
                                qn_ratio[d] = qn_ratio_f(lr_mag_in[d])
                                lr_LR[d] = lyra.LR(lr_fr[d],qn_ratio[d])
                            else:
                                qn_ratio[d] = 0
                                lr_LR[d] = 0
                            
                        #... temporarly memorize single ID's LRs:
                        bfactor[d] = lr_fr[d]
                        lr_single_outID.append(lr_LR[d])

                    elif bayes == True: # compute bayes factors and posterior probabilities
                        if (input_cat_type == 2) or (input_cat_type == 0):
                            if (output_cat_type == 2) or (output_cat_type == 0):
                                bpos = lyra.B_pos(sigma_output, sigma_input, Dra_tmp, Ddec_tmp) # astrometric bayes factor B_pos
                            else:
                                bpos = lyra.B_pos(sigma_output[i], sigma_input, Dra_tmp, Ddec_tmp)
                            if lr_mag_in[d] != nomag:
                                bphot = qn_ratio_f(lr_mag_in[d]) # photometric bayes factor -> LR
                            else:
                                qn_ratio[d] = 0
                                bphot = 1
                        
                        elif input_cat_type == 1:
                            if (output_cat_type == 2) or (output_cat_type == 0):
                                sigma_in_tmp = sigma_input[neigh_idx_lr[i]]
                                bpos = lyra.B_pos(sigma_output, sigma_in_tmp[j], Dra_tmp, Ddec_tmp)
                            else:
                                sigma_in_tmp = sigma_input[neigh_idx_lr[i]]
                                bpos = lyra.B_pos(sigma_output[i], sigma_in_tmp[j], Dra_tmp, Ddec_tmp)
                            if lr_mag_in[d] != nomag:
                                bphot = qn_ratio_f(lr_mag_in[d])
                            else:
                                qn_ratio[d] = 0
                                bphot = 1
                        
                        bfactor[d] = bpos * bphot * add_prior # bayes factor
                        ph = lyra.PH(dens1, dens2, Q_int) # prior
                        prior[d] = ph
                        phd = lyra.PHD(ph, bfactor[d]) # posterior
                        # instead of LR we have now the posterior:
                        lr_LR[d] = phd
                        lr_single_outID.append(lr_LR[d])

                    d+=1
            
            # from now on, everything is stored in lr_single_outID
            # compute Re with the LR-values of the single IDs matches:
            if len(lr_single_outID) != 0.:
                single_outID_LR = []
                for k in range(len(lr_single_outID)):
                    if bayes == False:
                        # rel.append(lyra.Re(np.sum(lr_single_outID), lr_single_outID[k], Q_int))
                        rel.append(lyra.Re(np.sum(lr_single_outID), lr_single_outID[k], bayes, Q_int))
                    else: # this is like p_i in N-WAY
                        # rel.append(lyra.Re(np.sum(lr_single_outID), lr_single_outID[k])) # completeness factor already included in the prior
                        rel.append(lyra.Re(np.sum(lr_single_outID), lr_single_outID[k], bayes))
                    single_outID_LR.append(lr_single_outID[k])
            
            else:
                lr_noID.append(ID_output[i])
                lr_noIDra.append(ra_output[i])
                lr_noIDdec.append(dec_output[i])
                if output_cat_type == 1:
                    sigma_output_tmp = sigma_output[i]
                else:
                    sigma_output_tmp = sigma_output
                lr_noIDerr.append(sigma_output_tmp)
                lr_noID_sec.append(sec_origin[i])  # this is fine to append tuples

            prog = 100*i/len(ID_output)
            print("Filling arrays: [%-20s] %d%%" % ('='*(int(prog/5)+1), (int(prog)+1)), end="\r")

        if len(lr_noID_sec)>0:
            assert all(isinstance(element, tuple) for element in lr_noID_sec), "Not all elements are tuples."
        if save_images == True:
            fig5, ax = plt.subplots()
            ax.scatter(lr_inRA, lr_inDEC, color='b', marker='o', alpha=0.5, label='Mathced input sources')
            ax.scatter(ra_output, dec_output, color='r', marker='.', label='All Output sources')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            ax.legend()
            ax.set_title('Matches ' + smag)
            ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
            plt.gca().invert_xaxis()
            fig5.savefig(path_images+'pos_match_results'+smag+add_name_output+'.png',bbox_inches="tight", dpi=200)

        print('\nPositional matches:', len(lr_LR))
        print('Total output sources not matched:',len(lr_noID))

        if save_output == True:
            lyra.print_output(path_output+filename_outinfo, 'Positional matches: '+ str(len(lr_LR)))
            lyra.print_output(path_output+filename_outinfo, 'Total output sources not matched: '+ str(len(lr_noID)))

    return globals()
