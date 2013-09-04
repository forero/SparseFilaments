	integer, parameter :: l = 8000000
	type (group), dimension(l) :: A
	real(dl), dimension(l) :: B
	integer, dimension(12) :: seed = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
	real(dl) :: random, x
	integer :: j, tmpi, orderB(l)
 
	    write (*,*) "Unsorted Values:"
	    call random_seed(put = seed)
	    do i = 1, l
	        call random_number(random)
	        A(i)%value = random
	        A(i)%order = i
	    end do
	 
	    call cpu_time(time1)
	    call QSort(A,l)
	    call cpu_time(time2)
	    write(*,*)  'Time used for quick sort : ', time2-time1
	    
	    write (*,*) "Sorted Values:"
	    do i = 1, 10
	        write (*,"((I8,1X,e14.7))") A(i*(l/10))
	    end do 
	    
	    !redo the same thing by using quick sort 2...
	    call random_seed(put = seed)
	    do i = 1, l
	        call random_number(random)
	        B(i) = random
	        orderB(i) = i
	    end do
	 
	    call cpu_time(time1)
	    call QSort2(B,orderB,l)
	    call cpu_time(time2)
	    write(*,*)  'Time used for quick sort2 : ', time2-time1
	    
	    write (*,*) "Sorted Values:"
	    do i = 1, 10
	        write (*,"(4(I8,1X,e14.7))") orderB(i*(l/10)), B(i*(l/10))
	    end do 
	    
	    !redo the same thing by using bubble sorting...
	    call random_seed(put = seed)
	    do i = 1, l
	        call random_number(random)
	        A(i)%value = random
	        A(i)%order = i
	    end do

	    call cpu_time(time1)
	    do j = 1, l
	    	do i = l, j+1, -1
	    		if(A(i-1)%value > A(i)%value) then
	    			x=A(i-1)%value; A(i-1)%value=A(i)%value; A(i)%value=x
	    			tmpi=A(i-1)%order; A(i-1)%order=A(i)%order; A(i)%order=tmpi
	    		endif
	    	enddo
	    enddo
	    call cpu_time(time2)
	    write(*,*)  'Time used for bubble sort : ', time2-time1
	    
	    
	    write (*,*) "Sorted Values:"
	    do i = 1, 10
	        write (*,"(4(I8,1X,e14.7))")  A(i*(l/10)-3:i*(l/10))
	    end do 
	    
	    stop
