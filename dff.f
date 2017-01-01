C FILE: dff.f
c      subroutine dff(SampleFrequency, spikes1, spikes2,
c     / Siteloc, max_ch1, max_ch2, purity, completeness, unit_match)

      subroutine dff(temp3, n_pixels,
     / Siteloc, max_ch1, max_ch2, purity, completeness, unit_match)

        integer spikes1(:,:), spikes2(:,:), array_temp(100000)
        real*4 dt, SampleFrequency, max_ptp_difference, max_distance
        real*4 purity(1000), completeness(1000)
        
        integer max_ch1(:), max_ch2(:), unit_match(1000,1)
        integer sort(1000,1000), sort_trans(1000,1000)
        integer comparematch(1000), k,j,p
        integer*2 Siteloc(:)  
        
Cf2py intent(in) SampleFrequency, spikes1, spikes2, max_ch1, max_ch2
Cf2py intent(in) Siteloc
Cf2py intent(out) purity, completeness, unit_match

        max_ch1 = max_ch1+1
        max_ch2 = max_ch2

        !dt=.5*SampleFrequency*1.0E-3    
        !max_distance=35.0           
        !max_distance=9999
        !max_ptp_difference = 0.50  
        !max_ptp_difference = 100
        
        !comparematch=999  !
        !sort = 0 !Keeps track of matching sort between sort1 and sort2



        for i in range(0,len(temp3),1):
            print "...averaging spike: ", i, " / ", len(temp3)
            #baseline = np.average(images_temp[temp3[i][0:len(temp3[i])/2]], axis=0) 
            baseline = np.average(images_temp[np.arange(temp3[i][0]-int(temp_img_rate)*3,temp3[i][0],1)], axis=0)  #Go back 3 sec prior to temp3[0] value...
            
            temp_img = images_temp[temp3[i]] #Remove avg of 1st half of images
            temp_frame = (temp_img - baseline)/baseline
            images_processed += temp_frame
        
        images_processed = images_processed/max(len(temp3),1)
            

        !loop over # ground-truth units
        do k = 1, size(spikes1,1)
            !write(*,*) "Cell: ", k
            !Load spikes into preexisting array; compute length
            do m = 1, size(spikes1(k,:))
                if (spikes1(k,m).eq.0) exit
                array_temp(m) = spikes1(k,m)
            enddo
            len_spikes1 = m-1
            !write(*,*) "Cell: ", k, " spks: ", len_spikes1

            !loop over # units in sorted data
            do j = 1, size(spikes2,1)
                x_dist = Siteloc(max_ch1(k)*2-1)-Siteloc(max_ch2(j)*2-1)
                y_dist = Siteloc(max_ch1(k)*2)-Siteloc(max_ch2(j)*2)
                unit_distance = sqrt(x_dist**2+y_dist**2)
                !write(*,*) x_dist, y_dist, unit_distance
                !Heuristic #1
                if (unit_distance.le.max_distance) then 
                  !write(*,*) "Compare to unit: ", j
                  do p = 1, len_spikes1
                     spike_temp = spikes2(j,p)  !spike p in sorted j
                     if (minval(abs(array_temp(1:len_spikes1)-
     &                     spike_temp)).le.dt) then
                         sort(k,j)=sort(k,j)+1
                     endif
                  enddo
                endif
            enddo

            !Save the bestmatching unit id
            !cell_match(k) = maxloc(sort(k,1:size(spikes2,1)))
            !write(*,*) "Matrix out: ", sort(k,1:size(spikes2,1))
c            write(*,*) "Max spike match: ", 
c     &            maxval(sort(k,1:size(spikes2,1))),
c     &            " unit: ", maxloc(sort(k,1:size(spikes2,1)))
        enddo

        !Transpose sort matrix to now look at units vs. cells
        sort_trans = transpose(sort)
        !sort_trans = sort
        !PURITY & COMPLETENESS COMPUTATION - Loop over sorted units
        do j = 1, size(spikes2,1)
            write(*,*) 
            write(*,*) "Unit: ", j
            !Compute exact size of spikes2 as it has trailing zeros
            do m = 1, size(spikes2(k,:))
                if (spikes2(j,m).eq.0) exit
            enddo
            len_spikes2 = m-1
            
            !Best matching cell; (j,:) vs (j,1) as maxloc is 2d array
            unit_match(j,:) = maxloc(sort_trans(j,1:size(spikes1,1)))
           
            !Purity: maxmatch_spks(unit)/tot_spks(unit)
            purity(j) = float(maxval(sort_trans(j,:)))/len_spikes2*100
                        
            !Compute exact size of cells(spikes1) of unit_match
            do m = 1, size(spikes1(unit_match(j,:),:))
                if (spikes1(unit_match(j,1),m).eq.0) exit
            enddo
            len_spikes1 = m-1
            
         write(*,*) "cell: ", unit_match(j,:), " #spks: ", len_spikes1
            !write(*,*) "size matching cell: ", len_spikes1
            !Completeness: maxmatch spks(unit)/tot_spks(matching cell)
            completeness(j) = float(maxval(sort_trans(j,:))) /
     &         float(len_spikes1)*100
            write(*,*) purity(j), completeness(j)
            if (completeness(j).ge.100) then
             write(*,*) "unit: ", j, " completeness: ", completeness(j)
            endif
            
        enddo
        !write(*,*) purity(1:size(spikes2,1))
        !write(*,*) completeness(1:size(spikes2,1))
        
c        return

      end subroutine compare_sort
C END FILE compare_sort.f



















