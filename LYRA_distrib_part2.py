
def LYRA_distrib_part2(args):
    globals().update(args)
    global save_output
    global save_images
    global path_output
    global filename_outinfo
    global r_in_tm
    global r_min_nm
    global r_max_nm
    global tree
    global nm_interp
    global qm_interp
    global r_out_tree
    global r_in_tree
    global r_lr
    global ra_input
    global dec_input
    global ra_output
    global dec_output
    global input_cat_type
    global ID_input
    global ID_output
    global mag_input
    global Q_int
    global Nout_int
    global area_output
    global area_input
    global area_int
    global needed_only_once
    global f
    global sigma_input
    global sigma_output
    global qn_ratio_f
    global exit_table_only
    global dens1
    global dens2
    global go_to_chull
    global polys_in
    global poly1
    global goto_pos_match
    global mag_failure


    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.path as mplPath
    import matplotlib
    matplotlib.use('Agg')
    from scipy import interpolate
    from scipy.spatial import KDTree
    import warnings
    import math
    import sys

    import LYRA_functions as lyra

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

    
    
    if needed_only_once == False:
        from LYRA_sigma_part1 import LYRA_sigma_part1
        globals().update(LYRA_sigma_part1(globals()))
        # outuput file info:
        if save_output == True:
            lyra.print_output(path_output+filename_outinfo, '\nSearching radii:\n'
                        + 'r_in_tn = ' + str(r_in_tm)
                        + ';  r_min_nm = ' +str(r_min_nm) + ';  r_max_nm = ' +str(r_max_nm)
                        + ';  r_lr = ' +str(r_lr))

        warnings.filterwarnings("ignore")

        ######################################################################### area calc ###
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
                            xint, yint = lyra.bartohistintersection(points_input_ra, points_input_dec, points_output_ra, points_output_dec)
                            
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

                if area_int > 0:
                    print(area_int)
                else:
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
                    lyra.print_output(path_output+filename_outinfo, '\n\n###### Error: No overlap between the catalogs. ######\n')

                exit_table_only = True
                from LYRA_pos_match_only import LYRA_pos_match_only
                LYRA_pos_match_only(globals())

        # source density
        dens1 = len(ra_input)/area_input
        dens2 = len(ra_output)/area_output
        print('Area primary and secondary catalog:', round(area_input,1), round(area_output,1), "deg2.")
        print('Overlapping area between catalogs:', round(area_int,3), 'deg2.')
        if save_output == True:
            lyra.print_output(path_output+filename_outinfo, 'Overlapping area between catalogs: '+str(round(area_int,2))+' deg2.')
        ############################################################################
        
        # KDTree implementation
        print("\nKDTree...")
        r_in_tree = np.column_stack((ra_input*3600, dec_input*3600))
        r_out_tree = np.column_stack((ra_output*3600, dec_output*3600))
        tree = KDTree(r_in_tree)
        print("... done.")

        needed_only_once = True


    ###################################
    ### Start working on magnitudes ###
    ###################################
    print("\n\n--- Work on mag_input:", smag, "---")
    if save_output == True:
        lyra.print_output(path_output+filename_outinfo, "\n--- Work on mag_input: " + smag + "---")

    ###################################################################### n(m) ########
    print('\nComputing n(m)...')
    neigh_idx_nmax = tree.query_ball_point(r_out_tree, r=r_max_nm, workers=-1)
    len_neigh_idx_nmax = tree.query_ball_point(r_out_tree, r=r_max_nm, workers=-1, return_length=True)
    len_neigh_idx_nmin = tree.query_ball_point(r_out_tree, r=r_min_nm, workers=-1, return_length=True)
    nm_sum = np.sum(len_neigh_idx_nmax)-np.sum(len_neigh_idx_nmin)

    if nm_sum < min_sample_size:
        print("\n### Warning: Not enough sources in the n(m) distribution, consider a different (r_min_nm, r_max_nm). Go to simple positional match. ###\n")
        if save_output == True:
            lyra.print_output(path_output+filename_outinfo, "### Warning: Not enough sources in the n(m) distribution, consider a different (r_min_nm, r_max_nm).  Go to simple positional match. ###")
        goto_pos_match = True
        mag_failure = True

    else:
        nn_mag_in = np.zeros(nm_sum)
        nn_in_ID = np.zeros(nm_sum)
        a=0
        for i in range(len(neigh_idx_nmax)):
            for j in range(len(neigh_idx_nmax[i])):
                rain_tmp = r_in_tree[neigh_idx_nmax[i],0]
                decin_tmp = r_in_tree[neigh_idx_nmax[i],1]
                Dra_nm = (r_out_tree[i,0] - rain_tmp[j])
                Ddec_nm = (r_out_tree[i,1] - decin_tmp[j])
                rr_nm = lyra.quadratic_sum(Dra_nm, Ddec_nm)
                if rr_nm > r_min_nm:
                    tmp_ID = ID_input[neigh_idx_nmax[i]]
                    nn_in_ID[a] = tmp_ID[j]
                    tmp_mag = mag_input[neigh_idx_nmax[i]]
                    nn_mag_in[a] = tmp_mag[j]
                    a=a+1

            prog = 100*i/len(neigh_idx_nmax)
            print("Filling arrays: [%-20s] %d%%" % ('='*(int(prog/5)+1), (int(prog)+1)), end="\r")

        unique_n, indices = np.unique(nn_in_ID, return_index=True)
        n_mag_in = np.zeros(len(indices))
        unique_n2 = np.zeros(len(indices))
        for i in range(len(indices)):
            n_mag_in[i] = nn_mag_in[indices[i]]
            unique_n2[i] = nn_in_ID[indices[i]]
        print('\n... done.')

        ###################################################################### q(m) ########
        print('\nComputing q(m)...')
        acircle = lyra.acircle(r_in_tm)
        aannu = lyra.aarea(r_min_nm, r_max_nm)

        #-------------------------------------------------------------- total(m):
        neigh_idx_totm = tree.query_ball_point(r_out_tree, r=r_in_tm, workers=-1)
        len_neigh_idx_totm = tree.query_ball_point(r_out_tree, r=r_in_tm, workers=-1, return_length=True)
        qm_sum = np.sum(len_neigh_idx_totm)

        if qm_sum < min_sample_size:
            print("\n### Warning: Not enough sources in the q(m) distribution, consider a different (r_in_tm). Go to simple positional match. ###\n")
            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, "### Warning: Not enough sources in the q(m) distribution, consider a different (r_in_tm).  Go to simple positional match. ###")
            goto_pos_match = True
            mag_failure = True
        
        else:

            tot_mag_in = np.zeros(qm_sum)
            tot_in_ID = np.zeros(qm_sum)
            a=c=f=0
            for i in range(len(neigh_idx_totm)): 
                b=0
                area_int = 0
                if go_to_chull == False:
                    for poly1 in polys_in:
                        if poly1.contains_points([[r_out_tree[i,0]/3600, r_out_tree[i,1]/3600]]) == True: #[ra_out, dec_out]
                            f+=1 # number of output sources inside the input area, i.e. inside also in the overlapping area
                        for j in range(len(neigh_idx_totm[i])):
                            rain_tmp = r_in_tree[neigh_idx_totm[i],0]
                            decin_tmp = r_in_tree[neigh_idx_totm[i],1]
                            Dra_tm = (r_out_tree[i,0] - rain_tmp[j])
                            Ddec_tm = (r_out_tree[i,1] - decin_tmp[j])
                            rr_tm = lyra.quadratic_sum(Dra_tm, Ddec_tm)
                            if rr_tm <= r_in_tm: # just a check, the if condition should be already verified
                                ID_tmp = ID_input[neigh_idx_totm[i]]
                                tot_in_ID[a] = ID_tmp[j]
                                magin_tmp = mag_input[neigh_idx_totm[i]]
                                tot_mag_in[a] = magin_tmp[j]
                                a+=1
                                b+=1
                    if b != 0:
                        c+=1 # at least one positional match between input and output sources
                else:
                    if poly1.contains_points([[r_out_tree[i,0]/3600, r_out_tree[i,1]/3600]]) == True: #[ra_out, dec_out]
                        f+=1 # number of output sources inside the input area, i.e. inside also in the overlapping area
                    for j in range(len(neigh_idx_totm[i])):
                        rain_tmp = r_in_tree[neigh_idx_totm[i],0]
                        decin_tmp = r_in_tree[neigh_idx_totm[i],1]
                        Dra_tm = (r_out_tree[i,0] - rain_tmp[j])
                        Ddec_tm = (r_out_tree[i,1] - decin_tmp[j])
                        rr_tm = lyra.quadratic_sum(Dra_tm, Ddec_tm)
                        if rr_tm <= r_in_tm: # just a check, the if condition should be already verified
                            ID_tmp = ID_input[neigh_idx_totm[i]]
                            tot_in_ID[a] = ID_tmp[j]
                            magin_tmp = mag_input[neigh_idx_totm[i]]
                            tot_mag_in[a] = magin_tmp[j]
                            a+=1
                            b+=1
                    if b != 0:
                        c+=1 # at least one positional match between input and output sources

                prog = 100*i/len(neigh_idx_totm)
                print("Filling arrays: [%-20s] %d%%" % ('='*(int(prog/5)+1), (int(prog)+1)), end="\r")

            unique_tot, indices = np.unique(tot_in_ID, return_index=True)
            total_mag_in = np.zeros(len(indices))
            unique_tot2 = np.zeros(len(indices))
            for i in range(len(indices)):
                total_mag_in[i] = tot_mag_in[indices[i]]
                unique_tot2[i] = tot_in_ID[indices[i]]
            print('\n... done.')

            # Nin = c
            # Nout = len(ID_output)
            #Q = 1 - (d/Nout) # because d/Nout = % of non-matched sources
            # estimation of Q parameter considering only sources in overlapping area
            # source_density_out = Nout/area_output
            # Nout_int = area_int * source_density_out
            #Q_int = 1 - (d/Nout_int)
            if (c > 0) and (f > 0):
                if compl_prior == True:
                    Q_int = c/f #(1-e/f) should give the same result
                else:
                    Q_int = 1
                print("Q_int", np.round(Q_int,6))
                # print(c,f)
            else:
                print('\n###### Error: No input sources within the chosen radius. ######\n')
                if save_output == True:
                    lyra.print_output(path_output+filename_outinfo, "###### Error: No input sources within the chosen radius. ######")

                exit_table_only = True
                from LYRA_pos_match_only import LYRA_pos_match_only
                LYRA_pos_match_only(globals())

            Nout_int = f
            print("\nTotal number of output sources:", len(ID_output))
            print("Total number of output sources inside the input area:", Nout_int)
            print("Completeness prior:", round(Q_int,6))

            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, "Total number of output sources: " + str(len(ID_output)))
                lyra.print_output(path_output+filename_outinfo, str(f) + " output sources inside the input area ("+str(round(f/len(ID_output),2))+"%)")
                lyra.print_output(path_output+filename_outinfo, 'Nin,matched = ' + str(int(c)) + ';  Nout,tot = ' + str(int(f)) + ';  Completeness prior = ' +str(round(Q_int,3)))

            print("\nComputing probability distributions...")
            bb=0
            for b in range(len(mag_input)):
                if mag_input[b] != nomag:
                    bb=bb+1
            tmp_mag_input = np.zeros(bb)
            bb=0
            for b in range(len(mag_input)):
                if mag_input[b] != nomag:
                    tmp_mag_input[bb] = mag_input[b]
                    bb=bb+1
            
            #rescaling n(m) area to tot(m) area:
            bins = np.linspace(min(tmp_mag_input),max(tmp_mag_input),distrib_bins)
            bar_width = bins[1]-bins[0]
            bin_height, bin_edges = np.histogram(n_mag_in, bins=bins)
            area_ratio = aannu/acircle
            bin_height = np.array(bin_height/area_ratio, dtype=float)
            bins_q = bins + (bins[1]-bins[0])/2
            bins_q = np.delete(bins_q,-1)

            # ... total(m):
            bin_hh, bin_e = np.histogram(total_mag_in, bins=bins, density=False)
            bin_h, bin_e = np.histogram(total_mag_in, bins=bins, density=True)
            # interpolate:
            Xplot_totm = np.linspace(min(tmp_mag_input), max(tmp_mag_input), len(bin_e)-1)
            funct_totm = interpolate.interp1d(bins_q, bin_h, bounds_error=False, fill_value='extrapolate')
            Yplot_totm = funct_totm(Xplot_totm)

            # ... background n(m):
            y_nm, x_nm = np.histogram(lyra.bartohist(bins_q, bin_height), density=True, bins=bins)
            bins_n = x_nm + abs(x_nm[0]-x_nm[1])/2
            Xplot_nm = Xplot_totm
            funct_nm = interpolate.interp1d(bins_n[:-1], y_nm, bounds_error=False, fill_value='extrapolate')
            Yplot_nm = funct_nm(Xplot_nm)
            Yplot_nm = np.where((np.isnan(Yplot_nm)==True), 0, Yplot_nm)

            # ... real(m):
            Xplot_realm = Xplot_totm
            Yplot_realm = Yplot_totm - Yplot_nm
            real_height = bin_hh - bin_height
            real_height = np.array(np.around(real_height,decimals=0), dtype=int)
            real_height = np.where(real_height<0, 0, real_height)
            real_hhist = np.around(bin_hh,0) - np.around(bin_height,0)
            real_xhist = lyra.bartohist(bins_q, real_hhist)
            y_realm, x_realm = np.histogram(real_xhist, density=True, bins=bins)
            funct_realm = interpolate.interp1d(bins_q, y_realm, bounds_error=False, fill_value='extrapolate')
            Yplot_realm = funct_realm(Xplot_realm)

            if save_images == True:
                fig4, ax = plt.subplots()
                ax.hist(total_mag_in, bins=bins, density=False, histtype='step', linestyle='--', linewidth=1.5, color='dodgerblue', label='total(m)')
                ax.hist(lyra.bartohist(bins_q, bin_height), color='mediumblue', density=False, histtype='step', linestyle='-.', linewidth=1.5, bins=bins, label='n(m)')
                YYreal, XXreal, patches = ax.hist(real_xhist, color='orange', density=False, histtype='stepfilled', alpha=0.5, bins=bins, label='real(m)')
                XXreal = XXreal + abs(XXreal[0]-XXreal[1])/2
                funct_realm_ = interpolate.interp1d(XXreal[:-1], YYreal, bounds_error=False, fill_value='extrapolate')
                func_to_plot = funct_realm_(np.linspace(min(Xplot_realm),max(Xplot_realm),100))
                ax.plot(np.linspace(min(Xplot_realm),max(Xplot_realm),100), np.where(func_to_plot<0, 0, func_to_plot), '-', color='k', label='real(m)')
                ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
                ax.set_ylabel('num. sources')
                ax.set_xlabel('Input magnitude/flux')
                ax.set_title('Counterparts distribution ' + smag)
                ax.legend(loc='upper left')
                fig4.savefig(path_images+'cparts_mag_distrib_'+smag+add_name_output+'.png',bbox_inches="tight", dpi=200)

            # q(m) = real(m)/somm(real(m))*Q normalize real(m) to the total source number
            #qm = (Yplot_realm/np.sum(y_realm))*Q_int
            qm = Yplot_realm*Q_int # normalizzation given by density=True in np.histogram.
            qm_interp = interpolate.interp1d(Xplot_realm, qm, bounds_error=False, fill_value='extrapolate')

            # same for n(m):
            nm = Yplot_nm
            nm_interp = interpolate.interp1d(Xplot_nm, nm, bounds_error=False, fill_value='extrapolate')

            #q(m)/n(m)
            q_ = qm_interp(Xplot_totm)
            n_ = nm_interp(Xplot_totm)
            n_ = np.where(n_<=0, 1, n_) # avoid division by zeros
            qn_ratio_f = interpolate.interp1d(Xplot_totm, q_/n_, bounds_error=False, fill_value='extrapolate')
            print("... done.")

    return globals()
