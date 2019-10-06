module compute_energy
  !--------------------------------------!
  !Input:
  ! pos, pos_ip0, pos_ip1, ip
  ! and the system parameters
  !Output: 
  ! EE, DeltaE 
  !--------------------------------------!
  implicit none

  save

!############coefficient in potential function#############!
  !
  !lj potential
  real*8,  private :: epsilon     !Energy unit epsilon in lj potential
  real*8,  private :: sigma       !Distance sigma in lj potential
  real*8,  private :: rc_lj       !Cut off radius of LJ potential
  real*8,  private :: rc_lj2      !rc_lj^2
  real*8,  private :: clx         !length of cell in x direction
  real*8,  private :: cly         !length of cell in y direction
  real*8,  private :: clz         !length of cell in z direction
  integer, private :: nclx        !number of cell in x direction
  integer, private :: ncly        !number of cell in y direction
  integer, private :: nclz        !number of cell in z direction 
!##########end coefficient in potential function###########!


!##########################arrays##########################!
  !
  !cell list in real space
  integer, allocatable, dimension( : ), private :: cell_list
  !
  !inverse cell list in real space
  integer, allocatable, dimension( : ), private :: inv_cell_list
  !
  ! head of chains, cell list
  integer, allocatable, dimension(:,:,:), private :: hoc    
  !
  ! head of chains, inverse cell list
  integer, allocatable, dimension(:,:,:), private :: inv_hoc
  !
  !neighbor cells of the center cell
  integer, allocatable, dimension(:,:,:), private :: cell_near_list 
!########################end arrays########################!


contains


subroutine initialize_energy_parameters
  !--------------------------------------!
  !Initial parameters are not inputted from file and compute
  !the total energy of the system.
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !Reference:
  !The computation of alpha, rc_real et al are refered to
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.304-306.
  !--------------------------------------!
  use global_variables
  implicit none
  !
  !read energy parameters from file
  call read_energy_parameters
  !
  !build lj_pair_list and lj_point
  call build_cell_list

end subroutine initialize_energy_parameters


subroutine read_energy_parameters
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  open(unit=100, file='energy_data.txt')
    read(100,*) epsilon
    read(100,*) sigma
    read(100,*) rc_lj
  close(100)

  nclx = int(Lx/rc_lj)     !cell numbers in x direction
  ncly = int(Ly/rc_lj)
  nclz = int(Lz/rc_lj)
  clx = Lx/nclx         !cell length    
  cly = Ly/ncly
  clz = Lz/nclz

  rc_lj2 = rc_lj*rc_lj

end subroutine read_energy_parameters


subroutine build_cell_list
  use global_variables
  implicit none
  integer :: i, j, k, l, m, n, p, q, r, x, y, z
  integer :: icelx, icely, icelz

  if (allocated(hoc)) deallocate(hoc)
  if (allocated(inv_hoc)) deallocate(inv_hoc)
  allocate(hoc(nclx,ncly,nclz))
  allocate(inv_hoc(nclx,ncly,nclz))
  hoc = 0
  inv_hoc = 0

  if (allocated(cell_list)) deallocate(cell_list)
  if (allocated(inv_cell_list)) deallocate(inv_cell_list)
  allocate(cell_list(NN))
  allocate(inv_cell_list(NN))
  cell_list = 0
  inv_cell_list = 0

  do i = 1, NN
    icelx = int((pos(i,1)-1)/clx) + 1
    icely = int((pos(i,2)-1)/cly) + 1
    icelz = int((pos(i,3)-1)/clz) + 1
    cell_list(i) = hoc(icelx,icely,icelz)
    hoc(icelx,icely,icelz) = i
  end do

  do i = NN, 1, -1
    icelx = int((pos(i,1)-1)/clx) + 1
    icely = int((pos(i,2)-1)/cly) + 1
    icelz = int((pos(i,3)-1)/clz) + 1
    inv_cell_list(i) = inv_hoc(icelx,icely,icelz)
    inv_hoc(icelx,icely,icelz) = i
  end do

  if(allocated(cell_near_list)) deallocate(cell_near_list)
  allocate(cell_near_list(nclx*ncly*nclz,28,3))
  cell_near_list = 0
  m = 0
  do i = 1, nclx
    do j = 1, ncly
      do k = 1, nclz
        m = m + 1
        n = 0
        do p = -1, 1
          do q = -1, 1
            do r = -1, 1
              x = i + p
              y = j + q
              z = k + r
              if (z>0 .and. z<=nclz) then
                n = n + 1
                if (x>nclx) then
                  x = x - nclx
                elseif (x<=0) then
                  x = x + nclx
                end if
                if (y>ncly) then
                  y = y - ncly
                elseif (y<=0) then
                  y = y + ncly
                end if
                cell_near_list(m,n,1) = x
                cell_near_list(m,n,2) = y
                cell_near_list(m,n,3) = z
              end if
            end do
          end do
        end do
        cell_near_list(m,28,1) = n
      end do
    end do
  end do

  open(100,file='./data/hoc_r.txt')
    do i = 1, nclx
      do j = 1, ncly
        do k = 1, nclz
         write(100,*) i,j,k,hoc(i,j,k),inv_hoc(i,j,k)
        end do
      end do
    end do
  close(100)

  open(100,file='./data/cell_list.txt')
    do i = 1, NN
      write(100,*) cell_list(i),inv_cell_list(i)
    end do
  close(100)

end subroutine build_cell_list


subroutine total_energy (EE)
  !--------------------------------------!
  !Compute total LJ potential energy,
  !including LJ energy of wall.
  !   
  !Input
  !   EE
  !Output
  !   EE
  !External Variables
  !   lj_point, lj_pair_list, pos, 
  !   epsilon, sigma, rc_lj, Lz
  !Routine Referenced:
  !1. rij_and_rr( rij, rr, i, j )
  !Reference:
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: EE
  integer :: i, j, k, icelx, icely, icelz, ncel
  real*8  :: rr, rij(3), inv_rr2, inv_rr6

  EE = 0

  do i = 1, NN
    icelx = int((pos(i,1)-1)/clx)+1
    icely = int((pos(i,2)-1)/cly)+1
    icelz = int((pos(i,3)-1)/clz)+1
    ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
    do j = 1, cell_near_list(ncel,28,1)
      icelx = cell_near_list(ncel,j,1)
      icely = cell_near_list(ncel,j,2)
      icelz = cell_near_list(ncel,j,3)
      k = hoc(icelx,icely,icelz)
      do while(k/=0)
        if (k/=i) then
          call rij_and_rr( rij, rr, i, k )
          if (rr<rc_lj2) then
            inv_rr2  = sigma*sigma/rr
            inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
            EE = EE + 4 * epsilon * ( inv_rr6 * inv_rr6 - inv_rr6 + 0.25D0) / 2
          end if
        end if
        k = cell_list(k)
      end do
    end do
  end do

end subroutine total_energy


subroutine Delta_Energy(DeltaE)
  !--------------------------------------!
  !Compute change of LJ potential Energy.
  !   
  !Input
  !   DeltaE
  !Output
  !   DeltaE
  !External Variables
  !   pos, lj_pair_list, lj_point
  !   pos_ip0, pos_ip1, ip
  !   Lx, Ly, Lz, sigma, epsilon, rc_lj
  !Routine Referenced:
  !
  !Reference:
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8  :: rij(3), rr, inv_rr2, inv_rr6
  real*8  :: EE1, EE2
  integer :: i,j,icelx,icely,icelz,ncel

  EE1 = 0
  icelx = int((pos_ip0(1)-1)/clx)+1
  icely = int((pos_ip0(2)-1)/cly)+1
  icelz = int((pos_ip0(3)-1)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc(icelx,icely,icelz) 
    do while (j/=0) 
      if (j/=ip) then
        call rij_and_rr(rij,rr,j,ip)
        if (rr<rc_lj2) then
          inv_rr2  = sigma*sigma/rr
          inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
          EE1 = EE1 + inv_rr6 * inv_rr6 - inv_rr6
        end if
      end if
      j = cell_list(j)
    end do
  end do
  EE1 = 4 * epsilon * EE1

  EE2 = 0
  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc(icelx,icely,icelz) 
    do while (j/=0) 
      if (j/=ip) then
        rij = pos_ip1 - pos(j,:)
        call periodic_condition1(rij)
        rr = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
        if (rr<rc_lj2) then
          inv_rr2  = sigma*sigma/rr
          inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
          EE2 = EE2 + inv_rr6 * inv_rr6 - inv_rr6
        end if
      end if
      j = cell_list(j)
    end do
  end do
  EE2 = 4 * epsilon * EE2

  DeltaE = EE2 - EE1

end subroutine Delta_Energy


subroutine update_cell_list
  !--------------------------------------!
  !Judge whether renew verlet list or not
  !   
  !Input
  !   EE
  !Output
  !   EE
  !External Variables
  !   Nq
  !Routine Referenced:
  !1.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: icelx, icely, icelz

  icelx = int((pos_ip0(1)-1)/clx)+1
  icely = int((pos_ip0(2)-1)/cly)+1
  icelz = int((pos_ip0(3)-1)/clz)+1
  call update_cell_list_delete(ip,icelx,icely,icelz)

  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
  call update_cell_list_add(ip,icelx,icely,icelz)

end subroutine update_cell_list


subroutine update_cell_list_delete(iq,icelx,icely,icelz)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer, intent(in) :: icelx, icely, icelz
  integer :: nti,bfi, j       ! next particle of ii, before particle of ii

  bfi = cell_list(iq)
  nti = inv_cell_list(iq)

  if ( bfi/=0 .and. nti/=0 ) then        !middle
    cell_list(nti) = bfi
    inv_cell_list(bfi) = nti
  elseif ( bfi==0 .and. nti/=0 ) then    !the first one
    cell_list(nti) = bfi
    inv_hoc(icelx,icely,icelz) = nti
  elseif ( bfi/=0 .and. nti==0 ) then    !the last one
    hoc(icelx,icely,icelz) = bfi
    inv_cell_list(bfi) = nti
  else                                   !only one
    hoc(icelx,icely,icelz) = nti
    inv_hoc(icelx,icely,icelz) = bfi
  end if

end subroutine update_cell_list_delete


subroutine update_cell_list_add(iq,icelx,icely,icelz)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer, intent(in) :: icelx, icely, icelz
  integer :: nti,bfi,j      ! next particle of ii, before particle of ii

  inv_cell_list(iq) = 0
  if ( inv_hoc(icelx,icely,icelz) /=0 ) then
    inv_cell_list( hoc(icelx,icely,icelz) ) = iq
  else
    inv_hoc(icelx,icely,icelz) = iq
  end if

  cell_list(iq) = hoc(icelx,icely,icelz)
  hoc(icelx,icely,icelz) = iq
  

end subroutine update_cell_list_add


subroutine compute_pressure (pressure)
  !----------------------------------------!
  !input:
  !  pos
  !output:
  !  pressure
  !External Variables:
  !  Ngl, Nml, Npe, NN,
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.52, Eq. (3.4.1).
  !----------------------------------------!
  use global_variables
  implicit none
  real*8, intent(out) :: pressure
  integer i,j,k
  real*8 :: rr, vir, inv_r2, inv_r6
  real*8, dimension(3) :: rij, fij

  vir = 0
  pressure = 0
  do i = 1, NN-1
    do j = i+1, NN
      call rij_and_rr(rij, rr, i, j)
      inv_r2 = sigma*sigma / rr
      inv_r6 = inv_r2*inv_r2*inv_r2
      fij = 48 * epsilon * inv_r2 * inv_r6 * (inv_r6-0.5) * rij
      vir = vir + dot_product(fij,rij)/3
    end do 
  end do
  pressure = rho / Beta + vir / (Lx*Ly*Lz) / NN
  write(*,*) pressure

end subroutine compute_pressure


end module compute_energy


















