C  Read in the complex refractive index for a material whose data is 
C     given in the file Material_refrac to resolution 1 or 2
C
      subroutine cmp_ref_index(Lambda_Res, Complex_refract1,
     *        Complex_refract2, Complex_refract3,
     *        Complex_refract4, Loss)
C
      implicit none
C
      integer*8 namelen, namelen2, Lambda_Res, i, Loss
      character Material_refrac*100, Material_refrac_b*100
      double precision Comp_ref1_tmp1(818), Comp_ref1_tmp2(818)
      double precision Comp_ref2_tmp1(1635), Comp_ref2_tmp2(1635)
      double precision rubbish
      complex*16 Complex_refract1(818), Complex_refract2(1635)
      complex*16 Complex_refract3(818), Complex_refract4(1635)
C
C
      open(3,file="Parameters/refrac_index.txt",status='old')
        read(3,'(a100)') Material_refrac
        read(3,'(a100)') Material_refrac_b
      close(3)
      namelen = len_trim(Material_refrac)
      namelen2 = len_trim(Material_refrac_b)
      if (Lambda_Res .eq. 1) then
        open(3,file="../Data/"//Material_refrac(1:namelen)//"_1.txt",
     *              status='unknown')
        open(4,file="../Data/"//Material_refrac_b(1:namelen2)//"_1.txt",
     *              status='unknown')
        if (Loss .eq. 1) then
          do i=1,818
            read(3,*) rubbish, Comp_ref1_tmp1(i), Comp_ref1_tmp2(i)
            Complex_refract1(i) = 
     *           dcmplx(Comp_ref1_tmp1(i), Comp_ref1_tmp2(i))
            read(4,*) rubbish, Comp_ref1_tmp1(i), Comp_ref1_tmp2(i)
            Complex_refract3(i) = 
     *           dcmplx(Comp_ref1_tmp1(i), 0.0d0)     ! only lossless substrate
          enddo
        else
          do i=1,818
            read(3,*) rubbish, Comp_ref1_tmp1(i), Comp_ref1_tmp2(i)
            Complex_refract1(i) = 
     *           dcmplx(Comp_ref1_tmp1(i), 0.0d0)
            read(4,*) rubbish, Comp_ref1_tmp1(i), Comp_ref1_tmp2(i)
            Complex_refract3(i) = 
     *           dcmplx(Comp_ref1_tmp1(i), 0.0d0)
          enddo
        endif
        close(3)
        close(4)
      elseif (Lambda_Res .eq. 2) then
        open(3,file="../Data/"//Material_refrac(1:namelen)//"_2.txt",
     *              status='unknown')
        open(4,file="../Data/"//Material_refrac_b(1:namelen2)//"_2.txt",
     *              status='unknown')
        if (Loss .eq. 1) then
          do i=1,1635
            read(3,*) rubbish, Comp_ref2_tmp1(i), Comp_ref2_tmp2(i)
            Complex_refract2(i) = 
     *           dcmplx(Comp_ref2_tmp1(i), Comp_ref2_tmp2(i))
            read(4,*) rubbish, Comp_ref2_tmp1(i), Comp_ref2_tmp2(i)
            Complex_refract4(i) = 
     *           dcmplx(Comp_ref2_tmp1(i), 0.0d0)     ! only lossless substrate
          enddo
        else
          do i=1,1635
            read(3,*) rubbish, Comp_ref2_tmp1(i), Comp_ref2_tmp2(i)
            Complex_refract2(i) = 
     *           dcmplx(Comp_ref2_tmp1(i), 0.0d0)
            read(4,*) rubbish, Comp_ref2_tmp1(i), Comp_ref2_tmp2(i)
            Complex_refract4(i) = 
     *           dcmplx(Comp_ref2_tmp1(i), 0.0d0)
          enddo
        endif
        close(3)
        close(4)
      endif
C
      return
      end 