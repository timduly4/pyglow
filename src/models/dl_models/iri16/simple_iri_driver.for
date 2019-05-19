
      program simple_iri_driver
          implicit none

          logical jf(50),print_output
          real alati,along,dhour,heibeg,heiend,heistp
          real oarr(100)
          real outf(20,1000)
          integer jmag,k
          integer mmdd,iyyyy

          !call read_ig_rz
          !call readapf107

          print_output = .true.

          do k=1,50
              jf(k) = .true.
          enddo

          jf(5)  = .false.
          jf(6)  = .false.
          jf(21) = .false.
          jf(22) = .false.
          jf(23) = .false.
          jf(29) = .false.
          jf(30) = .false.

          jmag = 0
          alati = 45.
          along = -80.
          iyyyy = 2002
          mmdd = -10
          dhour = 12.
          heibeg = 250.
          heiend = 251.
          heistp = 1.

          call read_ig_rz
          call readapf107

          call iri_sub(jf,jmag,alati,along,iyyyy,mmdd,dhour,
     &               heibeg,heiend,heistp,outf,oarr)


          if (print_output) then
              print *,'=============='
              print *,'  Outputs 1'
              print *,'=============='

              print *,'electron density, outf(1,*)     = ',outf(1,1)
              print *,'neutral temperature, outf(2,*)  = ',outf(2,1)
              print *,'ion temperature, outf(3,*)      = ',outf(3,1)
              print *,'electron temperature, outf(4,*) = ',outf(4,1)
              print *,'O+ ion density, outf(5,*)       = ',outf(5,1)
              print *,'H+ ion density, outf(6,*)       = ',outf(6,1)
              print *,'HE+ ion density, outf(7,*)      = ',outf(7,1)
              print *,'O2+ ion density, outf(8,*)      = ',outf(8,1)
              print *,'NO+ ion density, outf(9,*)      = ',outf(9,1)
          endif

          call iri_sub(jf,jmag,alati,along,iyyyy,mmdd,dhour,
     &               heibeg,heiend,heistp,outf,oarr)


          if (print_output) then
              print *,'=============='
              print *,'  Outputs 2'
              print *,'=============='

              print *,'electron density, outf(1,*)     = ',outf(1,1)
              print *,'neutral temperature, outf(2,*)  = ',outf(2,1)
              print *,'ion temperature, outf(3,*)      = ',outf(3,1)
              print *,'electron temperature, outf(4,*) = ',outf(4,1)
              print *,'O+ ion density, outf(5,*)       = ',outf(5,1)
              print *,'H+ ion density, outf(6,*)       = ',outf(6,1)
              print *,'HE+ ion density, outf(7,*)      = ',outf(7,1)
              print *,'O2+ ion density, outf(8,*)      = ',outf(8,1)
              print *,'NO+ ion density, outf(9,*)      = ',outf(9,1)
          endif

      end program simple_iri_driver
