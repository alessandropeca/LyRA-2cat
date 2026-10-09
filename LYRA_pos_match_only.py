
def LYRA_pos_match_only(args):
    globals().update(args)
    global save_output
    global path_output
    global filename_outinfo
    global ra_input
    global dec_input
    global ra_output
    global dec_output
    global ID_input
    global ID_output
    global sigma_input
    global sigma_output
    global r_lr
    global lr_noIDerr
    global sigma_match
    global exit_table_only
    global origin
    global lr_noID_sec
    global outID_names
    global inID_names
    
    from scipy.spatial import KDTree
    import numpy as np
    import pandas as pd
    import LYRA_functions as lyra
    import sys
    import warnings
    warnings.filterwarnings("ignore")

    if exit_table_only == False:
        print("\nPositional match only (no bayes):")
        if save_output == True:
            lyra.print_output(path_output+filename_outinfo, 'Positional match only (no bayes):')
    
    if needed_only_once == False:
        from LYRA_sigma_part1 import LYRA_sigma_part1
        globals().update(LYRA_sigma_part1(globals()))
    
    if (len(ra_input[ra_input != noxy]) > 0) and (len(ra_output[ra_output != noxy]) > 0):
        # KDTree implementation
        if exit_table_only == False:
            print("\nKDTree...")
        r_in_tree = np.column_stack((ra_input*3600, dec_input*3600))
        r_out_tree = np.column_stack((ra_output*3600, dec_output*3600))
        tree = KDTree(r_in_tree)
        if exit_table_only == False:
            print("... done.\n")

        neigh_idx = tree.query_ball_point(r_out_tree, r=r_lr, workers=-1)
        len_neigh_idx = tree.query_ball_point(r_out_tree, r=r_lr, workers=-1, return_length=True)
        neigh_sum = np.sum(len_neigh_idx)

        lr_inRA = np.zeros(neigh_sum)
        lr_inDEC = np.zeros(neigh_sum)
        lr_outRA = np.zeros(neigh_sum)
        lr_outDEC = np.zeros(neigh_sum)
        lr_inID = np.zeros(neigh_sum)
        lr_outID = np.zeros(neigh_sum)
        lr_r = np.zeros(neigh_sum)
        sigma_match = np.zeros(neigh_sum)
        # origin = np.empty(neigh_sum, dtype="<U9999")
        # origin = np.zeros(neigh_sum, dtype="<S9999")
        origin = [0] * neigh_sum
        s_flag = np.zeros(neigh_sum)
        lr_noID = []
        lr_noIDra = []
        lr_noIDdec = []
        lr_noIDerr = []
        lr_noID_sec = []
        
        d=0
        c_flag = []
        for i in range(len(neigh_idx)):
            lr_single_outID = []
            c_pos_match = []
            for j in range(len(neigh_idx[i])):
                rain_tmp = r_in_tree[neigh_idx[i],0]
                decin_tmp = r_in_tree[neigh_idx[i],1]
                Dra_tmp = (r_out_tree[i,0] - rain_tmp[j])
                Ddec_tmp = (r_out_tree[i,1] - decin_tmp[j])
                rr_lr = lyra.quadratic_sum(Dra_tmp, Ddec_tmp)
                if rr_lr <= r_lr:
                    ra_tmp = ra_input[neigh_idx[i]]
                    lr_inRA[d] = ra_tmp[j]
                    dec_tmp = dec_input[neigh_idx[i]]
                    lr_inDEC[d] = dec_tmp[j]
                    inID_tmp = ID_input[neigh_idx[i]]
                    lr_inID[d] = inID_tmp[j]
                    lr_outID[d] = ID_output[i]
                    lr_outRA[d] = ra_output[i]
                    lr_outDEC[d] = dec_output[i]
                    lr_r[d] = rr_lr
                    s_flag[d] = 1 # every source inside the matching radius is a possible match

                    if input_cat_type == 1:
                        sigma_input_tmp = sigma_input[j]
                    else:
                        sigma_input_tmp = sigma_input
                    if output_cat_type == 1:
                        sigma_output_tmp = sigma_output[i]
                    else:
                        sigma_output_tmp = sigma_output
                    sigma_match[d] = min(sigma_output_tmp, sigma_input_tmp)
                    # pri_origin_tmp = pri_origin[neigh_idx[i]]
                    pri_origin_tmp = [pri_origin[i_] for i_ in neigh_idx[i]]

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

                    # if isinstance(pri_origin_tmp[j], str) == False:
                    #     pri_origin_tmp[j] = str(pri_origin_tmp[j])
                    # if isinstance(sec_origin[i], str) == False:
                    #     sec_origin[i] = str(sec_origin[i])
                    # underscore_str = "_"
                    # origin[d] = pri_origin_tmp[j] + underscore_str.encode() + sec_origin[i]
                    
                    # pri_origin_tmp_str = pri_origin_tmp[j]
                    # if isinstance(pri_origin_tmp_str, str) == False:
                    #     pri_origin_tmp_str = str(pri_origin_tmp_str)
                    # sec_origin_str = sec_origin[i]
                    # if isinstance(sec_origin_str, str) == False:
                    #     sec_origin_str = str(sec_origin_str)
                    # underscore_str = "_"
                    # pri_origin_tmp_str.replace('lyra_', '').replace('.csv', '')
                    # sec_origin_str.replace('lyra_', '').replace('.csv', '')
                    # origin[d] = (pri_origin_tmp_str + underscore_str + sec_origin_str).encode()
                
                    lr_single_outID.append(lr_r[d])
                    c_pos_match.append(rr_lr) # local matching radii
                    d=d+1
            
            if len(lr_single_outID) == 0.:
                lr_noID.append(ID_output[i])
                lr_noIDra.append(ra_output[i])
                lr_noIDdec.append(dec_output[i])
                if output_cat_type == 1:
                    sigma_output_tmp = sigma_output[i]
                else:
                    sigma_output_tmp = sigma_output
                lr_noIDerr.append(sigma_output_tmp)
                lr_noID_sec.append(sec_origin[i])

            else:
                for ff in range(len(c_pos_match)):
                    if c_pos_match[ff] == min(c_pos_match):
                        c_flag.append(1)
                    else:
                        c_flag.append(0)
            
            prog = 100*i/len(ID_output)
            #if exit_table_only == False:
                #print("Filling array: [%-20s] %d%%" % ('='*(int(prog/5)+1), (int(prog)+1)), end="\r")

        if len(lr_noID_sec)>0:
            assert all(isinstance(element, tuple) for element in lr_noID_sec), "Not all elements are tuples."

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
        sigma_match = np.append(sigma_match,lr_noIDerr)
        origin.extend(lr_noID_sec)
        lr_inRA = np.append(lr_inRA, np.full((len(lr_noID),),-99))
        lr_inDEC = np.append(lr_inDEC, np.full((len(lr_noID),),-99))
        lr_r = np.append(lr_r, np.full((len(lr_noID),),-99))
        s_flag = np.append(c_flag, np.full((len(lr_noID),),-99))
        c_flag = np.append(c_flag, np.full((len(lr_noID),),-99))
        c_flag = np.where(c_flag!=1, 0, c_flag)

        
        data_to_write = {"OutID": outID_names,
                        # "OutID": lr_outID,
                        "RAout": lr_outRA,
                        "DECout": lr_outDEC,
                        "InID": inID_names,
                        # "InID": lr_inID,
                        "RAin": lr_inRA,
                        "DECin": lr_inDEC,
                        "match_err": sigma_match,
                        "origin": origin,
                        "r": lr_r,
                        "S_flag": s_flag,
                        "C_flag": c_flag
                        }
        
        df = pd.DataFrame(data_to_write)
        df = df.sort_values(by=["OutID"])

        print('\n\nComputing flags...')
        df_ok = df[df["InID"]!= -99]
        s_outID = np.unique(df_ok['OutID'])
        index=0
        for subid in s_outID:
            index+=1
            df_sub = df[df['OutID']==subid]
            cf_tmp = list(df_sub['C_flag'])
            b_index = df_sub['C_flag'].index.tolist()
            if cf_tmp.count(1) > 1:
                list_r_tmp = df_sub['r'].tolist()
                best_r_idx = list_r_tmp.index(max(list_r_tmp))
                best_index = b_index[best_r_idx]
                try:
                    if len(best_index) > 1: # it means we have two entries with the same max value
                        best_index = best_index[0] # so both the sources are ok, take one
                except:
                    best_index = best_index
                
                indeces = df_sub.index.tolist()
                for idx in indeces:
                    if idx == best_index:
                        df['C_flag'][idx] = 1
                    else:
                        df['C_flag'][idx] = 0
            prog = 100*index/len(s_outID)
            # print("[%-20s] %d%%" % ('='*(int(prog/5)+1), (int(prog)+1)), end="\r")
        # print('\n\nComputing flags...')
        # df_ok = df[df["InID"] != -99]
        # s_outID = np.unique(df_ok['OutID'])

        # for index, subid in enumerate(s_outID, start=1):
        #     df_sub = df[df['OutID'] == subid]
        #     cf_tmp = df_sub['C_flag'].tolist()
        #     b_index = df_sub['C_flag'].idxmax()
            
        #     if cf_tmp.count(1) > 1:
        #         list_r_tmp = df_sub['r'].tolist()
        #         best_r_idx = list_r_tmp.index(max(list_r_tmp))
        #         best_index = b_index[best_r_idx]
                
        #         if isinstance(best_index, list) and len(best_index) > 1:
        #             best_index = best_index[0]
                    
        #         indeces = df_sub.index.tolist()
        #         df['C_flag'].iloc[indeces] = np.where(indeces == best_index, 1, 0)
            
        #     prog = 100 * index / len(s_outID)
        #     # print("[%-20s] %d%%" % ('='*(int(prog/5)+1), (int(prog)+1)), end="\r")


        # print('\n... done.')
        print('\nCleaning multiple in ID counterpars...')
        df_sc = df.sort_values(by=["OutID"])
        df_sc_ok = df_sc[df_sc["InID"]!= -99]
        s_inID = np.unique(df_sc_ok['InID'])
        check_multiple_inID = df_sc_ok[df_sc_ok['C_flag'] == 1]['InID'].duplicated().any() # returns True if more than one inID with C_flag == 1
        if check_multiple_inID == True:
            index=0
            for subid in s_inID:
                index+=1
                df_sub = df_sc[df_sc['InID']==subid]
                df_sub = df_sub[df_sub['C_flag']==1] ####
                cf_tmp = list(df_sub['C_flag'])
                b_index = df_sub['C_flag'].index.tolist()
                
                if cf_tmp.count(1) > 1:
                    list_r_tmp = df_sub['r'].tolist()
                    r_val3 = min(list_r_tmp)
                    r_best_idx = df_sub[df_sub['r'] == r_val3].index
                    if len(r_best_idx)>1: # same minimum -> both are ok
                        r_best_idx = r_best_idx[0]
                    for idx in b_index:
                        if idx == r_best_idx:
                            outid_tmp = df_sc['OutID'][idx]
                            def_sub2 = df_sc[(df_sc['OutID'] == outid_tmp) & (df_sc['C_flag'] == 1)]
                            if len(def_sub2) <= 1:
                                df_sc['C_flag'][idx] = 1
                            else:
                                df_sc['C_flag'][idx] = 0
                            df_sc['C_flag'][idx] = 1
                        else:
                            df_sc['C_flag'][idx] = 0
                
                prog = 100*index/len(s_inID)
                # print("[%-20s] %d%%" % ('='*(int(prog/5)), (int(prog))), end="\r")
        print('\n... done.')
        # print('\nCleaning multiple in ID counterparts...')
        # df_sc = df.sort_values(by=["OutID"])
        # df_sc_ok = df_sc[df_sc["InID"] != -99]
        # s_inID = np.unique(df_sc_ok['InID'])
        # check_multiple_inID = df_sc_ok[df_sc_ok['C_flag'] == 1]['InID'].duplicated().any()

        # if check_multiple_inID:
        #     for index, subid in enumerate(s_inID, start=1):
        #         df_sub = df_sc[(df_sc['InID'] == subid) & (df_sc['C_flag'] == 1)]
        #         cf_tmp = df_sub['C_flag'].tolist()
        #         b_index = df_sub.index.tolist()
                
        #         if cf_tmp.count(1) > 1:
        #             list_r_tmp = df_sub['r'].tolist()
        #             r_val3 = min(list_r_tmp)
        #             r_best_idx = df_sub[df_sub['r'] == r_val3].index[0]
                    
        #             for idx in b_index:
        #                 if idx == r_best_idx:
        #                     outid_tmp = df_sc['OutID'].iloc[idx]
        #                     def_sub2 = df_sc[(df_sc['OutID'] == outid_tmp) & (df_sc['C_flag'] == 1)]
        #                     if len(def_sub2) <= 1:
        #                         df_sc['C_flag'].iloc[idx] = 1
        #                     else:
        #                         df_sc['C_flag'].iloc[idx] = 0
        #                     df_sc['C_flag'].iloc[idx] = 1
        #                 else:
        #                     df_sc['C_flag'].iloc[idx] = 0
                
        #         prog = 100 * index / len(s_inID)
        #         # print("[%-20s] %d%%" % ('='*(int(prog/5)), (int(prog))), end="\r")

        
        # set origin of the matched catalog
        file_pri_tw_a = file_primary.split("/")[-1]
        if file_pri_tw_a[0] == 'd':
            file_pri_tw_ = file_pri_tw_a.split(".csv")[0]
            file_pri_tw = (int(file_pri_tw_.split("d")[-1]),)
        else:
            file_pri_tw_b = file_pri_tw_a.split("lyra_")[-1]
            file_pri_tw = (int(file_pri_tw_b.split(".csv")[0]),)
        
        file_sec_tw_a = file_secondary.split("/")[-1]
        if file_sec_tw_a[0] == 'd':
            file_sec_tw_ = file_sec_tw_a.split(".csv")[0]
            file_sec_tw = (int(file_sec_tw_.split("d")[-1]),)
        else:
            file_sec_tw_b = file_sec_tw_a.split("lyra_")[-1]
            file_sec_tw = (int(file_sec_tw_b.split(".csv")[0]),)
        
        # df["origin"] = np.where(df_sc["InID"] != -99, file_pri_tw + file_sec_tw, file_sec_tw)
        assert all(isinstance(element, tuple) for element in df_sc["origin"].to_list()), "Not all elements are tuples."
        df_sc["origin"] = [file_pri_tw + file_sec_tw if inid != -99 else file_sec_tw for inid in df_sc["InID"]]
        assert all(isinstance(element, tuple) for element in df_sc["origin"].to_list()), "Not all elements are tuples."
        
        df = df_sc.sort_values(by=["OutID"])
        print('\n... done.')
    
        best_matches = df['C_flag'].tolist()
        if exit_table_only == False:
            print("\n\nTotal matches:", d)
            print("Best matches:", best_matches.count(1))
            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, "Total matches: " + str(d))
                lyra.print_output(path_output+filename_outinfo, "Best matches: " + str(best_matches.count(1)))

        # write LYRA outputs
        out_file_name = path_output+filename_LR+".csv"
        df.to_csv(out_file_name, sep=',', index=False, header=True)
        if exit_table_only == True:
            print('Exit.')
            sys.exit()
    else:
        print('\n###### Error: No input/output sources with good coordinates. ######\n')
        if save_output == True:
            lyra.print_output(path_output+filename_outinfo, '\n###### No input/output sources with good coordinates. ######')
        
        if (len(ra_output[ra_input != noxy]) > 0):
            
            exit_err = np.sqrt(ra_output_err**2 + dec_output_err**2)
            lr_outID = lr_outID.astype(np.int32)
            outID_names = ID_output_names[lr_outID]
            data_to_write = {"OutID": outID_names,
                        # "OutID": ID_output,
                        "RAout": ra_output,
                        "DECout": dec_output,
                        "match_err": exit_err,
                        "r": 0,
                        "S_flag": 0,
                        "C_flag": 0
                        }
        
            df = pd.DataFrame(data_to_write)
            df = df.sort_values(by=["OutID"])
            # write LYRA outputs
            out_file_name = path_output+filename_LR+".csv"
            df.to_csv(out_file_name, sep=',', index=False, header=True)
            print('Exit.')
            sys.exit() # to exit from LYRA

        else:
            print('Exit.')
            sys.exit() # to exit from LYRA