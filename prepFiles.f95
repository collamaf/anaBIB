program dump_processor

    implicit double precision (a-h,o-z)
    implicit integer (i-n)
    character :: flag1,flag2
    INTEGER :: stat
    character(10) :: elem,elem2
    character(80) :: outputFileName0
    integer, dimension(1000) :: ifin

    open(30,file='part_new_sigma_25_nocut')
    open(40,file='info_new_sigma_25_nocut')

    acheck=0
    nsimtot=0
    do j=1,74
      write(outputFileName0, '(a,I0,a)') 'sigma250/bkg_IP2_bpipe_sigma_250m_collim',j,'_DUMP'
      open(20,file= outputFileName0,status='unknown')

      i=1
      do while (i>0)
        read(20,*,IOSTAT=stat) flag1
        IF(IS_IOSTAT_END(stat)) then
          ifin(j)=i-1
          write(40,*) ifin(j)
          backspace 20
          backspace 20
          read(20,*) flag2, flag2,flag2,flag2, nsim
          i=-20
        end if
        i=i+1
      end do
      close(20)
      nsimtot=nsimtot+nsim
    end do



    do j=1,74
      write(outputFileName0, '(a,I0,a)') 'sigma250/bkg_IP2_bpipe_sigma_250m_collim',j,'_DUMP'
      open(20,file= outputFileName0,status='unknown')
      do ii=1,ifin(j)
        read(20,*) flag1
        if (ichar(flag1) .ne.35) then
          backspace 20
          !read(20,fmt='( 1x, I5, 5(1x,E12.5), 3(1x,F11.3) , E12.5, 1X, A8, 3(1x,F11.3), I2, 1X, A8,  I4)') &
          !&  npartsec, a1, a2, a3, a4, a5, a6,a7,a8,a9,elem,a10,a11,a12,na,elem2,na2
          read(20,*)  npartsec, a1, a2, a3, a4, a5, a6,a7,a8,a9,elem,a10,a11,a12,na,elem2,na2
          if ( abs(npartsec)==211 .or. abs(npartsec)==321) npartsec=-2212
          a5=a5*2e12/nsimtot
          !if (a9 < 30) then
          !if (a4 < 0) then
          if(a12<2500) then
              if(npartsec==2112) then
                !if(a1>0.0001) then
                  write(30,fmt='(i8,9es16.6,a16,3es16.6,i8,a16,i8)') npartsec,a1, a2, a3, a4, a5,a6,a7,&
                  a8,a9,elem,a10,a11,a12,na,elem2,na2
                  acheck=acheck+1
              !  end if
              else if (abs(npartsec)==2212) then
              !  if(a1>0.001) then
                  write(30,fmt='(i8,9es16.6,a16,3es16.6,i8,a16,i8)') npartsec,a1, a2, a3, a4, a5,a6,a7,&
                  a8,a9,elem,a10,a11,a12,na,elem2,na2
                  acheck=acheck+1
              !  end if
              else if (abs(npartsec)==22) then
              !  if(a1>0.0002) then
                  write(30,fmt='(i8,9es16.6,a16,3es16.6,i8,a16,i8)') npartsec,a1, a2, a3, a4, a5,a6,a7,&
                  a8,a9,elem,a10,a11,a12,na,elem2,na2
                  acheck=acheck+1
              !  end if
              else if (abs(npartsec)==11) then
              !  if(a1>0.0002) then
                  write(30,fmt='(i8,9es16.6,a16,3es16.6,i8,a16,i8)') npartsec,a1, a2, a3, a4, a5,a6,a7,&
                  a8,a9,elem,a10,a11,a12,na,elem2,na2
                  acheck=acheck+1
              !  end if
              else
                write(30,fmt='(i8,9es16.6,a16,3es16.6,i8,a16,i8)') npartsec,a1, a2, a3, a4, a5,a6,a7,a8,a9,elem,&
                a10,a11,a12,na,elem2,na2
                acheck=acheck+1
              end if
            !end if
          end if
        end if
      end do
      close(20)
    end do

    print*,acheck
    close(30)
    close(40)
    stop
end
