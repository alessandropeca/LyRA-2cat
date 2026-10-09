
def LYRA_pos_match_only_bayes(args):
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
    global r_lr
    global sigma_input
    global sigma_output
    global input_cat_type
    global output_cat_type
    global lr_noIDerr
    global sigma_match
    global tree
    global exit_table_only
    global rel_threshold
    global idx_to_change
    global origin
    global lr_noID_sec
    global mag_failure
    global outID_names
    global inID_names
    
    from scipy.spatial import KDTree
    import numpy as np
    import pandas as pd
    import matplotlib.path as mplPath
    import matplotlib.pyplot as plt
    import LYRA_functions as lyra
    import sys
    import warnings
    warnings.filterwarnings("ignore")

    def clockwiseangle_and_distance(point): # to order points clockwise
        origin = xy_centroid
        refvec = [0,1]
        # Vector between point and the origin: v = p - o
        vector = [point[0]-origin[0], point[1]-origin[1]]
        # Length of vector: ||v||
        lenvector = math.hypot(vector[0], vector[1])
        # If length is zero there is no angle
        if lenvector == 0:
            return -math.pi, 0
        # Normalize vector: v/||v||
        normalized = [vector[0]/lenvector, vector[1]/lenvector]
        dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
        diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
        angle = math.atan2(diffprod, dotprod)
        # Negative angles represent counter-clockwise angles so we need to subtract them from 2*pi (360 degrees)
        if angle < 0:
            return 2*math.pi+angle, lenvector
        # I return first the angle because that's the primary sorting criterium
        # but if two vectors have the same angle then the shorter distance should come first.
        return angle, lenvector

    print("\nPositional match only:")
    if save_output == True:
        lyra.print_output(path_output+filename_outinfo, '\nPositional match only:')

    if mag_failure == False:
        from LYRA_sigma_part1 import LYRA_sigma_part1
        globals().update(LYRA_sigma_part1(globals()))
    
    # Area calc:
    print('Area check...')
    go_to_chull = False
    if convhull_area == False: # clustering (time consuming but more reliable)
        try:
            polys_in, area_input = lyra.find_pgons(ra_input, dec_input) # polygons here are mpl Paths
            polys_out, area_output = lyra.find_pgons(ra_output, dec_output)
            area_int = 0
            for poly1 in polys_in:
                points_input_ra, points_input_dec = poly1.vertices[:, 0], poly1.vertices[:, 1]

                for poly2 in polys_out:
                    if (poly1.intersects_path(poly2, filled=True)) == True: # with filled=True, it returns True even if one polygon is inside the other
                        
                        points_output_ra, points_output_dec = poly2.vertices[:, 0], poly2.vertices[:, 1]
                        # check for area overlap
                        xint, yint = lyra.intersection(points_input_ra, points_input_dec, points_output_ra, points_output_dec)
                        
                        if len(xint) == 0: # one polygon is inside the other
                            if area_input > area_output:
                                area_int = area_output
                            else:
                                area_int = area_input
                        
                        else: # polygons intersection
                            # find intersection polygon
                            xy_intersect = np.column_stack((xint, yint))
                            
                            for v in range(len(points_input_ra)):
                                if poly2.contains_point((points_input_ra[v],points_input_dec[v])) == True:
                                    xy_intersect = np.append(xy_intersect, [[points_input_ra[v],points_input_dec[v]]], axis=0)
                            for v in range(len(points_output_ra)):
                                if poly1.contains_point((points_output_ra[v],points_output_dec[v])) == True:
                                    xy_intersect = np.append(xy_intersect, [[points_output_ra[v],points_output_dec[v]]], axis=0)
                            xy_centroid = [sum(xy_intersect[:,0])/len(xy_intersect[:,0]), sum(xy_intersect[:,1])/len(xy_intersect[:,1])]
                            
                            xy_intersect = np.array(sorted(xy_intersect, key=clockwiseangle_and_distance))
                            xy_intersect = np.append(xy_intersect, [xy_intersect[0]], axis=0) #close the polygon
                            area_int = area_int + lyra.PolyArea(xy_intersect[:,0], xy_intersect[:,1])

            if area_int <= 0:
                print('###### No overlap between the catalogs. ######\n')
                area_int = 0.
                if save_output == True:
                    lyra.print_output(path_output+filename_outinfo, '\n###### Error: No overlap between the catalogs. ######')

                exit_table_only = True
                from LYRA_pos_match_only import LYRA_pos_match_only
                LYRA_pos_match_only(globals())
        except:
            error = sys.exc_info()[0]
            print("Warning: an error occured in clustering procedure, go to convex hull estimation.")
            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, '\n###### Error: '+str(error)+' ######')
            go_to_chull = True

    else:
        go_to_chull = True

    if go_to_chull == True:
        points_input_ra, points_input_dec = lyra.convex_hull(ra_input,dec_input)
        points_output_ra, points_output_dec = lyra.convex_hull(ra_output,dec_output)
        # compute areas
        area_input = lyra.PolyArea(points_input_ra, points_input_dec)
        area_output = lyra.PolyArea(points_output_ra, points_output_dec)
        #Polygon creation poly1
        poly1 = mplPath.Path(np.stack((points_input_ra, points_input_dec), axis=1))
        poly2 = mplPath.Path(np.stack((points_output_ra, points_output_dec), axis=1))

        # check areas overlap
        if (poly1.intersects_path(poly2, filled=True)) == True: # with filled=True, it returns True even if one polygon is inside the other
            # check for area overlap
            xint, yint = lyra.intersection(points_input_ra, points_input_dec, points_output_ra, points_output_dec)
            
            if len(xint) == 0: # one polygon is inside the other
                if area_input > area_output:
                    area_int = area_output
                else:
                    area_int = area_input
            
            else: # polygons intersection
                # find intersection polygon
                xy_intersect = np.column_stack((xint, yint))
                
                for v in range(len(points_input_ra)):
                    if poly2.contains_point((points_input_ra[v],points_input_dec[v])) == True:
                        xy_intersect = np.append(xy_intersect, [[points_input_ra[v],points_input_dec[v]]], axis=0)
                for v in range(len(points_output_ra)):
                    if poly1.contains_point((points_output_ra[v],points_output_dec[v])) == True:
                        xy_intersect = np.append(xy_intersect, [[points_output_ra[v],points_output_dec[v]]], axis=0)
                xy_centroid = [sum(xy_intersect[:,0])/len(xy_intersect[:,0]), sum(xy_intersect[:,1])/len(xy_intersect[:,1])]
                xy_intersect = np.array(sorted(xy_intersect, key=clockwiseangle_and_distance))
                xy_intersect = np.append(xy_intersect, [xy_intersect[0]], axis=0) #close the polygon
                area_int = lyra.PolyArea(xy_intersect[:,0], xy_intersect[:,1])
        else:
            print('\n###### Error: No overlap between the catalogs. ######\n')
            area_int = 0.
            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, '\n###### Error: No overlap between the catalogs. ######')

            exit_table_only = True
            from LYRA_pos_match_only import LYRA_pos_match_only
            LYRA_pos_match_only(globals())
    
    # source density
    dens1 = len(ra_input)/area_input
    dens2 = len(ra_output)/area_output
    print('Area primary and secondary catalog:', round(area_input,1), round(area_output,1), "deg2.")
    print('Overlapping area between catalogs:', round(area_int,2), 'deg2.')
    if save_output == True:
        lyra.print_output(path_output+filename_outinfo, 'Overlapping area between catalogs: '+str(round(area_int,2))+' deg2.')
    print("... done.")

    # Matching:
    if (len(ra_input[ra_input != noxy]) > min_sample_size) and (len(ra_output[ra_output != noxy]) > min_sample_size):
        # KDTree implementation
        print("\nKDTree...")
        r_in_tree = np.column_stack((ra_input*3600, dec_input*3600))
        r_out_tree = np.column_stack((ra_output*3600, dec_output*3600))
        tree = KDTree(r_in_tree)
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
        lr_LR = np.zeros(neigh_sum)
        bfactor = np.zeros(neigh_sum)
        prior = np.zeros(neigh_sum)
        sigma_match = np.zeros(neigh_sum)
        # origin = np.empty(neigh_sum, dtype="<U9999")
        # origin = np.zeros(neigh_sum, dtype="<S9999")
        origin = [0] * neigh_sum
        lr_noID = []
        lr_noIDra = []
        lr_noIDdec = []
        lr_noIDerr = []
        lr_noID_sec = []

        c=f=0
        for i in range(len(neigh_idx)):
            b=0
            area_int = 0
            if convhull_area == False:
                for poly1 in polys_in:
                    if poly1.contains_points([[r_out_tree[i,0]/3600, r_out_tree[i,1]/3600]]) == True: #[ra_out, dec_out]
                        f+=1 # number of output sources inside the overlapping area
                    for j in range(len(neigh_idx[i])):
                        rain_tmp = r_in_tree[neigh_idx[i],0]
                        decin_tmp = r_in_tree[neigh_idx[i],1]
                        Dra_tmp = (r_out_tree[i,0] - rain_tmp[j])
                        Ddec_tmp = (r_out_tree[i,1] - decin_tmp[j])
                        rr_lr = lyra.quadratic_sum(Dra_tmp, Ddec_tmp)
                        if rr_lr <= r_lr:
                            b+=1
                if b != 0:
                    c+=1 # at least one positional match between input and output sources
                
                prog = 100*i/len(ID_output)
                #print("Compute distributions: [%-20s] %d%%" % ('='*(int(prog/5)+1), (int(prog)+1)), end="\r")
            else:
                if poly1.contains_points([[r_out_tree[i,0]/3600, r_out_tree[i,1]/3600]]) == True: #[ra_out, dec_out]
                    f+=1 # number of output sources inside the overlapping area
                for j in range(len(neigh_idx[i])):
                    rain_tmp = r_in_tree[neigh_idx[i],0]
                    decin_tmp = r_in_tree[neigh_idx[i],1]
                    Dra_tmp = (r_out_tree[i,0] - rain_tmp[j])
                    Ddec_tmp = (r_out_tree[i,1] - decin_tmp[j])
                    rr_lr = lyra.quadratic_sum(Dra_tmp, Ddec_tmp)
                    if rr_lr <= r_lr:
                        b+=1
                if b != 0:
                    c+=1 # at least one positional match between input and output sources
                
                prog = 100*i/len(ID_output)
                #print("Compute distributions: [%-20s] %d%%" % ('='*(int(prog/5)+1), (int(prog)+1)), end="\r")
        
        # Nout = len(ID_output)
        # estimation of Q parameter considering only sources in overlapping area
        # source_density_out = Nout/area_output
        # Nout_int = area_int * source_density_out
        if (c > 0) and (f > 0):
            if compl_prior == True:
                Q_int = c/f
            else:
                Q_int = 1
            
            print("\nTotal number of output sources:", len(ID_output))
            print("Total number of output sources inside the input area:", f)
            print("Completeness prior:", np.round(Q_int,6))
            
            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, "Total number of output sources: " + str(len(ID_output)))
                lyra.print_output(path_output+filename_outinfo, str(f) + " output sources inside the input area ("+str(round(f/len(ID_output),2))+"%)")
                lyra.print_output(path_output+filename_outinfo, 'Nin,matched = ' + str(int(c)) + ';  Nout,tot = ' + str(int(f)) + ';  Completeness prior = ' +str(round(Q_int,6)))
        else:
            print('\n###### Error: No input sources within the chosen radius. ######\n')
            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, '\n###### Error: No input sources within the chosen radius. ######')

            exit_table_only = True
            from LYRA_pos_match_only import LYRA_pos_match_only
            LYRA_pos_match_only(globals())
            
        Nout_int = f
        d=0
        rel = []
        #pany = []
        for i in range(len(neigh_idx)):
            lr_single_outID = []
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
                    #... temporarly memorize single ID's LRs:
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

                    # bayesian statistics:
                    if (input_cat_type == 2) or (input_cat_type == 0):
                        if (output_cat_type == 2) or (output_cat_type == 0):
                            bpos = lyra.B_pos(sigma_output, sigma_input, Dra_tmp, Ddec_tmp) # astrometric bayes factor B_pos 
                        else:
                            bpos = lyra.B_pos(sigma_output[i], sigma_input, Dra_tmp, Ddec_tmp)
                    
                    elif input_cat_type == 1:
                        if (output_cat_type == 2) or (output_cat_type == 0):
                            sigma_in_tmp = sigma_input[neigh_idx[i]]
                            bpos = lyra.B_pos(sigma_output, sigma_in_tmp[j], Dra_tmp, Ddec_tmp)
                        else:
                            sigma_in_tmp = sigma_input[neigh_idx[i]]
                            bpos = lyra.B_pos(sigma_output[i], sigma_in_tmp[j], Dra_tmp, Ddec_tmp)
                    
                    bfactor[d] = bpos * add_prior
                    ph = lyra.PH(dens1, dens2, Q_int)
                    prior[d] = ph
                    phd = lyra.PHD(ph, bfactor[d])
                    lr_LR[d] = phd
                    lr_single_outID.append(lr_LR[d])
                    d+=1

            # compute reliabilities and pany
            if len(lr_single_outID) != 0.:
                single_outID_LR = []
                single_outID_rel = []
                #pany_tmp = lyra.Pany(lr_single_outID)
                for k in range(len(lr_single_outID)):
                    rel_tmp = lyra.Re(np.sum(lr_single_outID), lr_single_outID[k], bayes) # no need Q_int here, since it is aleady implemented in ph
                    rel.append(rel_tmp)
                    single_outID_LR.append(lr_single_outID[k])
                    single_outID_rel.append(rel_tmp)
                    #pany.append(pany_tmp)

            else:
                lr_noID.append(ID_output[i])
                lr_noIDra.append(ra_output[i])
                lr_noIDdec.append(dec_output[i])
                if output_cat_type == 1:
                    sigma_output_tmp = sigma_output[i]
                else:
                    sigma_output_tmp = sigma_output
                lr_noIDerr.append(sigma_output_tmp)
                lr_noID_sec.append(sec_origin[i]) # this is fine to append tuples
            
                
            prog = 100*i/len(ID_output)
            #print("Filling array: [%-20s] %d%%" % ('='*(int(prog/5)+1), (int(prog)+1)), end="\r")
                
        if len(lr_noID_sec)>0:
            assert all(isinstance(element, tuple) for element in lr_noID_sec), "Not all elements are tuples."
        print("\n\nComputing thresholds...")
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

        print("... done.")
        LR_threshold = LRth[RCmax_index]

        print('\nPth: ' + str('{:.4g}'.format(LRth[RCmax_index])))
        print('C[Pth]: ' + str('{:.4g}'.format(C[RCmax_index])))
        print('R[Pth]: ' + str('{:.4g}'.format(R[RCmax_index])))

        if save_images == True:
            fig6, ax = plt.subplots()
            ax.plot(LRth, RC, color='r', label='(R+C)/2', zorder=3)
            ax.plot(LRth, R, color='mediumblue', ls='--', label='R', zorder=3)
            ax.plot(LRth, C, color='mediumblue', ls=':', label='C', zorder=3)
            ax.axvline(color='r', lw=3, ls='-.', alpha=0.4, x=LRth[RCmax_index])
            ax.set_xlabel('$LR_{th}$')
            ax.grid(ls=':', color='grey', alpha=0.3, zorder=0)
            ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
            ax.set_title('Completeness Bayes')
            # ax.set_xscale('log')
            ax.legend()
            fig6.savefig(path_images+'completeness_bayes.png',bbox_inches="tight", dpi=200)

        
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
        origin.extend(lr_noID_sec) # to concatenate two lists of tuples
        lr_inRA = np.append(lr_inRA, np.full((len(lr_noID),),-99))
        lr_inDEC = np.append(lr_inDEC, np.full((len(lr_noID),),-99))
        lr_r = np.append(lr_r, np.full((len(lr_noID),),-99))
        bfactor = np.append(bfactor, np.full((len(lr_noID),),-99))
        lr_LR = np.append(lr_LR, np.full((len(lr_noID),),-99))
        prior = np.append(prior, np.full((len(lr_noID),),-99))
        rel = np.append(rel, np.full((len(lr_noID),),-99))

        data_to_write = {"OutID": outID_names,
                        # "OutID": lr_outID,
                        "RAout": lr_outRA,
                        "DECout": lr_outDEC,
                        "InID": inID_names,
                        # "InID": lr_inID,
                        "RAin": lr_inRA,
                        "DECin": lr_inDEC,
                        "r": lr_r,
                        "match_err": sigma_match,
                        "origin": origin,
                        "Bfact": np.log10(bfactor),
                        "Prior": prior,
                        "Post": lr_LR,
                        #"Pany": pany,
                        "Prel": rel,
                        "F_flag": -99,
                        "S_flag": 0,
                        "C_flag": 0
                        }

        print()
        print(len(outID_names), len(origin))
        df = pd.DataFrame(data_to_write)
        df = df.sort_values(by=["OutID"])
        

        print('\n\nComputing flags...')
        # F_flag (Flavour_flag): 1 = best match with LR>=LRth, 2 = secondary match with LR>=LRth and rel>=rel_threshold, 0 = other cases LR<LRth, -99 = no match.
        # C_flag: 1 = best match, 0 = low prob. or no match.
        df_ok = df[df["InID"]!= -99]
        s_outID = np.unique(df_ok['OutID'])
        index=0
        for subid in s_outID:
            index+=1
            df_sub = df[df['OutID']==subid] # df_sub has same indeces of df
            lr_tmp = df_sub['Post']
            prel_tmp = df_sub['Prel']

            if len(df_sub[df_sub.Post >= LR_threshold]) > 0:
                idx_list = df_sub.index.to_list()
                for idx in idx_list:
                    max_lr = max(lr_tmp)
                    if (lr_tmp[idx] == max_lr) and (lr_tmp[idx] >= LR_threshold):
                        df['F_flag'][idx] = 1
                        df['C_flag'][idx] = 1
                    else:
                        if (lr_tmp[idx] >= LR_threshold) and (prel_tmp[idx] >= rel_threshold):
                            df['F_flag'][idx] = 2
                        else:
                            df['F_flag'][idx] = 0
                        df['C_flag'][idx] = 0
            
            prog = 100*index/len(s_outID)
            # print("[%-20s] %d%%" % ('='*(int(prog/5)), (int(prog))), end="\r")
        
        # print('\n\nComputing flags...')
        # df_ok = df[df["InID"] != -99]
        # s_outID = np.unique(df_ok['OutID'])
        # for index, subid in enumerate(s_outID, start=1):
        #     df_sub = df[df['OutID'] == subid]  # df_sub has same indices of df
        #     lr_tmp = df_sub['Post']
        #     prel_tmp = df_sub['Prel']

        #     if len(df_sub[df_sub.Post >= LR_threshold]) > 0:
        #         idx_list = df_sub.index.to_list()
        #         max_lr = max(lr_tmp)
        #         df['F_flag'].iloc[idx_list] = (lr_tmp == max_lr) & (lr_tmp >= LR_threshold)
        #         df['C_flag'].iloc[idx_list] = df['F_flag'].iloc[idx_list]

        #         df['F_flag'].iloc[idx_list] = np.where((lr_tmp >= LR_threshold) & (prel_tmp >= rel_threshold), 2, df['F_flag'].iloc[idx_list])
        #         df['C_flag'].iloc[idx_list] = np.where((lr_tmp >= LR_threshold) & (prel_tmp >= rel_threshold), 0, df['C_flag'].iloc[idx_list])
            
        #     prog = 100 * index / len(s_outID)
        #     # print("[%-20s] %d%%" % ('='*(int(prog/5)), (int(prog))), end="\r")
        
        # b = len(df[df['C_flag']==1])
        # print(str(b) + '/' + str(f) + '  (' + str(np.round((100*b/f),2)) + '%) output sources with a good match (P >= Pth) in the overalapping area.\n')

        # if save_output == True:
        #     lyra.print_output(path_output+filename_outinfo, '(R+C)/2 max: ' + str(round(RC[RCmax_index],2)) + '\n'
        #                 + 'Pth: ' + str('{:.3g}'.format(LRth[RCmax_index])) + '\n'
        #                 + 'C[Pth]: ' + str('{:.3g}'.format(C[RCmax_index])) + '\n'
        #                 + 'R[Pth]: ' + str('{:.3g}'.format(R[RCmax_index])) + '\n'
        #                 + str(b) + '/' + str(f) + '  (' + str(np.round((100*b/f),2)) + '%) output sources with a good match (P >= Pth) in the overalapping area.')
        
        print('... done.')
        print('\nCleaning multiple inID counterpars...')
        df_sc = df.sort_values(by=["OutID"])
        df_sc_ok = df_sc[df_sc["InID"]!= -99]
        s_inID = np.unique(df_sc_ok['InID'])

        check_multiple_inID = df_sc_ok[df_sc_ok['C_flag'] == 1]['InID'].duplicated().any() # returns True if more than one inID with C_flag == 1
        if check_multiple_inID == True:
            index=0
            for subid in s_inID:
                index+=1
                df_sub = df_sc[df_sc['InID']==subid]
                df_sub = df_sub[df_sub['C_flag']==1]
                b_index = df_sub['C_flag'].index.tolist()
                if len(df_sub) > 1: # if more than one inID has C_flag == 1, then find the one with the best LR/posterior
                    list_LR_tmp = df_sub['Post'].tolist()
                    LR_best_idx = df_sub[df_sub['Post'] == max(list_LR_tmp)].index
                    try:
                        if len(LR_best_idx) > 1: # it means we have two entries with the same max value
                            LR_best_idx = LR_best_idx[0] # so both the sources are ok, take one
                    except:
                        LR_best_idx = LR_best_idx
                    for idx in b_index:
                        if idx == LR_best_idx:
                            df_sc['C_flag'][idx] = 1 # should be already == 1
                        else:
                            df_sc['C_flag'][idx] = 0 # we need to check if now there are reliable inID for the corresponding outID which are not already C_flag == 1
                            
                            df_out_tmp = df_sc[(df_sc['OutID'] == df_sc['OutID'][idx]) & (df_sc['OutID'] != -99)] # what is the outID corresponding to that idx? df_sc['OutID'][idx]
                            clist = df_out_tmp[df_out_tmp['C_flag'] == 1]
                            if len(clist) < 1: # there is not a C_flag == 1 for this outID. check if F_flag == 2
                                df_out_f2_idx = df_out_tmp[df_out_tmp['F_flag'] == 2].index.to_list()
                                id_list_tmp = df_sc[df_sc['C_flag']==1]
                                id_list = id_list_tmp['InID'].tolist()
                                if (len(df_out_f2_idx) == 1) and (df_sc['InID'][df_out_f2_idx[0]] not in id_list): # this new inID must not be already takes
                                    df_sc['C_flag'][df_out_f2_idx[0]] = 1
                                elif (len(df_out_f2_idx) > 1) and (df_sc['InID'][df_out_f2_idx[0]] not in id_list): # then check the max reliability
                                    df_out_f2 = df_out_tmp[df_out_tmp['F_flag'] == 2]
                                    rel2 = df_out_f2['Prel']
                                    for idx2 in df_out_f2_idx:
                                        if rel2[idx2] == max(rel2):
                                            df_sc['C_flag'][idx2] = 1

                    prog = 100*index/len(s_inID)
                    #print("[%-20s] %d%%" % ('='*(int(prog/5)), (int(prog))), end="\r")
        # print('\nCleaning multiple inID counterparts...')
        # df_sc = df.sort_values(by=["OutID"])
        # df_sc_ok = df_sc[df_sc["InID"] != -99]
        # s_inID = np.unique(df_sc_ok['InID'])

        # check_multiple_inID = df_sc_ok[df_sc_ok['C_flag'] == 1]['InID'].duplicated().any()

        # if check_multiple_inID:
        #     for index, subid in enumerate(s_inID, start=1):
        #         df_sub = df_sc[(df_sc['InID'] == subid) & (df_sc['C_flag'] == 1)]
        #         b_index = df_sub.index.tolist()

        #         if len(df_sub) > 1:
        #             list_LR_tmp = df_sub['Post'].tolist()
        #             LR_best_idx = df_sub[df_sub['Post'] == max(list_LR_tmp)].index[0]

        #             for idx in b_index:
        #                 if idx == LR_best_idx:
        #                     df_sc.at[idx, 'C_flag'] = 1
        #                 else:
        #                     df_sc.at[idx, 'C_flag'] = 0

        #                     df_out_tmp = df_sc[(df_sc['OutID'] == df_sc.at[idx, 'OutID']) & (df_sc['OutID'] != -99)]
        #                     clist = df_out_tmp[df_out_tmp['C_flag'] == 1]

        #                     if len(clist) < 1:
        #                         df_out_f2_idx = df_out_tmp[df_out_tmp['F_flag'] == 2].index.tolist()

        #                         if len(df_out_f2_idx) == 1:
        #                             df_sc.at[df_out_f2_idx[0], 'C_flag'] = 1
        #                         elif len(df_out_f2_idx) > 1:
        #                             df_out_f2 = df_out_tmp[df_out_tmp['F_flag'] == 2]
        #                             rel2 = df_out_f2['Prel']

        #                             for idx2 in df_out_f2_idx:
        #                                 if rel2[idx2] == max(rel2):
        #                                     df_sc.at[idx2, 'C_flag'] = 1

        #         prog = 100 * index / len(s_inID)
        #         #print("[%-20s] %d%%" % ('=' * (int(prog/5)), (int(prog))), end="\r")


        print('\n... done.')
        # S_flag == 1 if F_flag == 1.
        s_flag = np.where(np.array(df_sc['F_flag'])==1, 1, 0)
        df_sc['S_flag'] = s_flag
        column_to_move = df_sc.pop('S_flag')
        df_sc["S_flag"] = column_to_move
        column_to_move = df_sc.pop('C_flag')
        df_sc["C_flag"] = column_to_move

        # set origin of the matched catalog
        # file_pri_tw_a = file_primary.split("/")[-1]
        # if file_pri_tw_a[0] == 'd':
        #     file_pri_tw_ = file_pri_tw_a.split(".csv")[0]
        #     file_pri_tw = (int(file_pri_tw_.split("d")[-1]),)
        # else:
        #     file_pri_tw_b = file_pri_tw_a.split("lyra_")[-1]
        #     file_pri_tw = (int(file_pri_tw_b.split(".csv")[0]),)
        
        # file_sec_tw_a = file_secondary.split("/")[-1]
        # if file_sec_tw_a[0] == 'd':
        #     file_sec_tw_ = file_sec_tw_a.split(".csv")[0]
        #     file_sec_tw = (int(file_sec_tw_.split("d")[-1]),)
        # else:
        #     file_sec_tw_b = file_sec_tw_a.split("lyra_")[-1]
        #     file_sec_tw = (int(file_sec_tw_b.split(".csv")[0]),)

        # df_sc["origin"] = np.where(df_sc["InID"] != -99, file_pri_tw + file_sec_tw, file_sec_tw)

        assert all(isinstance(element, tuple) for element in df_sc["origin"].to_list()), "Not all elements are tuples."
        
        # df_sc["origin"] = [file_pri_tw + file_sec_tw if inid != -99 else file_sec_tw for inid in df_sc["InID"]]
        
        assert all(isinstance(element, tuple) for element in df_sc["origin"].to_list()), "Not all elements are tuples."
        
        df = df_sc.sort_values(by=["OutID"])
        best_matches = df['C_flag'].tolist()
        print("\n\nPositional matches:", d)
        print("Best matches:", best_matches.count(1))
        if save_output == True:
            lyra.print_output(path_output+filename_outinfo, "-----------------------------------\nPositional matches: " + str(d))
            lyra.print_output(path_output+filename_outinfo, "Best matches: " + str(best_matches.count(1)))

        # write LYRA outputs
        out_file_name = path_output+filename_LR+".csv"
        df.to_csv(out_file_name, sep=',', index=False, header=True)
        
        # print('\nCheck duplicates:')
        # out_dupl, in_dupl = lyra.check_duplicates(path_output+filename_LR+".csv")
        # print('Out duplicates:', out_dupl)
        # print('In duplicates:', in_dupl)
        # print()
        
    else:
        print('\n###### Error: No input/output sources with good coordinates. ######\n')
        if save_output == True:
            lyra.print_output(path_output+filename_outinfo, '\n###### No input/output sources with good coordinates. ######')

        exit_table_only = True
        from LYRA_pos_match_only import LYRA_pos_match_only
        LYRA_pos_match_only(globals())