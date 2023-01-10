!Author: H. Paul Keeler, 2020.
!Website: hpaulkeeler.com
!Repository: github.com/hpaulkeeler/posts
!For more details, see post:
!hpaulkeeler.com/simulating-poisson-random-variables-in-fortran/

program test_poisson
    implicit none
 
   !lambda is the Poisson parameter (that is, its mean)
   real,parameter :: lambda = 4.7 
   !number of variables
   integer :: numbSim = 1000 
 
   real  :: meanPoisson !(sample) mean
   real :: varPoisson !(sample) variance
 
   !declare functions
   integer funPoissonSingle 
   real funUniformSingle
   
   !START Collect statistists on Poisson variables
   !initialize statistics
   integer :: numbPoissonTemp
   real :: sumPoisson = 0
   real :: sumPoissonSquared = 0
   
   !loop through for each random variable
   integer :: i !count (loop) variable
   do i = 1, numbSim
      !generate a single poisson variable
       numbPoissonTemp = funPoissonSingle(lambda)
 
       !total sum of variables
       sumPoisson = sumPoisson+numbPoissonTemp
       !total sum of squared variables
       sumPoissonSquared =sumPoissonSquared+ numbPoissonTemp**2
 
       if (i <= 5) then
          !print the first 5 numbers
          print *, "One of the Poisson variables has the value ", numbPoissonTemp
       end if
       
     end do
     
     !calculate statistics
     meanPoisson = sumPoisson / numbSim 
     varPoisson = sumPoissonSquared / numbSim - meanPoisson**2 
 
     !print statistics
     print *, "The average of the Poisson variables is ", meanPoisson
     print *, "The variance of the Poisson variables is", varPoisson
     print *, "For Poisson random variables, the mean and variance will agree more and more as the number of simulations increases."
 
     !END Collect statistists on Poisson variables
 end program test_poisson
 
 !START Function definitions
 !Uniform function -- returns a standard uniform random variable
 function funUniformSingle() result(randUni)
     implicit none
     real randUni
     call random_seed
     call random_number(randUni) 
 end function
 
 !Poisson function -- returns a single Poisson random variable
 function funPoissonSingle(lambda) result(randPoisson)
     implicit none
     real, intent(in) :: lambda !input
     real :: exp_lambda !constant for terminating loop
     real :: randUni !uniform variable 
     real :: prodUni !product of uniform variables
     integer :: randPoisson !Poisson variable
 
     !declare functions
     real funUniformSingle
     
     exp_lambda= exp(-lambda) !calculate constant
 
     !initialize variables
     randPoisson = -1
     prodUni = 1
     do while (prodUni > exp_lambda)       
          randUni = funUniformSingle() !generate uniform variable
          prodUni = prodUni * randUni !update product
          randPoisson = randPoisson + 1 !increase Poisson variable
     end do  
 end function
 !END  Function definitions