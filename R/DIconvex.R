
DIconvex<-function(x, lower, upper,increasing=FALSE, epsim=0, epsic=0,visual=TRUE ){


  #sanity checks fot x
  stopifnot(!is.null(x))
  stopifnot(is.finite(x))
  stopifnot(is.numeric(x))


  if (anyDuplicated(x) != 0){
    stop("No duplicated data is allowed in the vector x")
  }

  if(is.unsorted(x) == TRUE){
    stop("The data in the vector x has to be ranked in ascending order")
  }

  #sanity checks for lower
  stopifnot(!is.null(lower))
  stopifnot(is.finite(lower))
  stopifnot(is.numeric(lower))

  #sanity checks for upper
  stopifnot(!is.null(upper))
  stopifnot(is.finite(upper))
  stopifnot(is.numeric(upper))
  stopifnot(all(lower<upper))

  #other sanity checks
  stopifnot(length(lower)== length(upper) & length(x) == length(upper))


  #increasing function
  DIconvex_I<-function(x, lower, upper,epsim, epsic, visual){



    data = cbind(x, lower, upper)
    puts = as.data.frame(data)

    G<-2*(max(upper)-min(lower))


    p_xl<-nrow(puts)                #sets the number of puts
    p_fl<-p_xl                      #the number of put flags
    p_zl<-p_xl                      #flag i times put price i



    # x1,...,xlen are the prices, we turn on off certain options if they do not satisfy increasing or convexity constraint
    # introduce  auxiliary variables
    # zpp1 := p2 fp1,...,pmi fpmi-1, zcp1:=c
    p_vl<-p_xl+p_fl+p_zl            #Number of decision variables

    p_mylp<-make.lp(2*p_xl # correspondence x_i and z_i
                    +p_xl-1 # monotonicity
                    +p_xl-2 # convexity
                    +2* p_xl # the lower and upper bounds for the y which we need for the new formulation
                    ,p_vl)        #Creates a new linear program model object.Specifies the number of
    #constraints (nrow) and the number of decision variables (ncol).

    #mylp<-make.lp(4*zpl,vl)
    #row.add.mode(mylp,state="on")

    set.objfn(p_mylp,obj=rep(1,p_fl),indices=(p_xl+1):(p_xl+p_fl))           #Sets the objective function. Specifies: the linear program model object,
    #the coefficients of the objective function (in this case a vector with as many 1s as flags
    #and 0s elsewhere).

    lp.control(p_mylp,sense="max")                                           #Specifies that it is a maximization problem.

    set.type(p_mylp,(p_xl+1):(p_xl+p_fl),type="binary")                      #Sets the type of the decision variables to "binary" (= "integer" with an upper
    #bound of 1 and lower bound of 0). In our case the decision variable can assume the value of 1 or 0.

    #Sets the type of the constraints, which are all less or equal ("<=")
    set.constr.type(p_mylp, c(rep(1,2*p_xl+p_xl-1+p_xl-2+2* p_xl)))

    p_iecount<-1
    p_rhs<-NULL

    #Here we compile the matrix of our linear program model object. We consider x binary, L lower bound (lower), U upper bound (upper).
    # additional definitions to not make too many mistakes with the indices
    i_x     <- 0
    i_fl    <- p_xl
    i_zl    <- p_xl+p_fl

    #P_zpl:= f_i*x_(i+1). Flag i times put price i+1
    # z<=y-L(1-f)
    # z>=y-U(1-f)

    # z-y-Lf<=-L
    # -z+Uf+y<=U


    # f1 x2
    # first two constraints with flagbar, second two with flags
    for (i in 1:p_xl){

      set.row(p_mylp,p_iecount, xt=c(1,-1,G),indices=c(i_zl+i, i_x+i, i_fl+i))    #Third constraint corresponding to y+Ux-z
      p_iecount <- p_iecount + 1

      set.row(p_mylp, p_iecount,xt=c(-1,1,G), indices=c(i_zl+i,i_x+i,i_fl+i))  #Fourth constraint corresponding to z-Lx-y
      p_iecount <- p_iecount + 1

      p_rhs<-c(p_rhs,c(G,G))                                     #The right side of the constraints corresponding to the vector [0,0, U, -L]
    }

    # are options increasing in x? x_i - x_(i+1) <= 0.
    # f1 (x2-x1)

    for (i in 2:(p_xl)){
      set.row(p_mylp, p_iecount,xt=c(-1,1), indices=c(i_zl+i,i_zl+i-1))
      p_iecount <- p_iecount + 1

      p_rhs<-c(p_rhs,-epsim)

    }

    # are options convex in x?
    # O(K2)<=alpha O(K1)+(1-alpha) O(K3)
    # we do it with flags here
    # f2*(O(K2)-alpha O(k1)-(1-alpha)O(K3))<=0

    for (i in 2:(p_xl-1)){
      alpha<-(puts$x[i]-puts$x[i+1])/(puts$x[i-1]-puts$x[i+1])
      set.row(p_mylp, p_iecount,xt=c(1,-alpha,alpha-1),
              indices=c(i_zl+i,i_zl+i-1,i_zl+i+1))
      p_iecount <- p_iecount + 1

      p_rhs<-c(p_rhs,-epsic)

    }



    #lower and upper bounds for x_i
    for (i in 1:p_xl){

      set.row(p_mylp,p_iecount, xt=c(-1),indices=c( i_x+i))    #Third constraint corresponding to y+Ux-z
      p_iecount <- p_iecount + 1

      set.row(p_mylp, p_iecount,xt=c(1), indices=c( i_x+i))  #Fourth constraint corresponding to z-Lx-y
      p_iecount <- p_iecount + 1

      p_rhs<-c(p_rhs,c(-puts$lower[i],puts$upper[i]))                                     #The right side of the constraints corresponding to the vector [0,0, U, -L]
    }


    set.rhs(p_mylp,p_rhs)                          #sets the right hand side of the linear program model object


    # lp.control(p_mylp,bb.depthlimit=0)


    p_sol<-solve(p_mylp)


    p_flags<-as.logical(get.variables(p_mylp)[(p_xl+1):(p_xl+p_fl)])
    # print(p_flags)


    solvec_dec<-get.variables(p_mylp)[(p_xl+p_fl+1):(p_xl+p_fl+p_xl)]

    if(visual){
      matplot(puts$x, cbind(puts$lower,solvec_dec,puts$upper),type='l')
    }

    return(list(  p_flags, solvec_dec[p_flags],p_sol))

  }

  #decreasing function

  DIconvex_D<-function(x, lower, upper,epsim, epsic,visual){




    data = cbind(x, lower, upper)
    puts = as.data.frame(data)

    G<-2*(max(upper)-min(lower))


    p_xl<-nrow(puts)                #sets the number of puts
    p_fl<-p_xl                      #the number of put flags
    p_zl<-p_xl                      #flag i times put price i



    # x1,...,xlen are the prices, we turn on off certain options if they do not satisfy increasing or convexity constraint
    # introduce  auxiliary variables
    # zpp1 := p2 fp1,...,pmi fpmi-1, zcp1:=c
    p_vl<-p_xl+p_fl+p_zl            #Number of decision variables

    p_mylp<-make.lp(2*p_xl # correspondence x_i and z_i
                    +p_xl-1 # monotonicity
                    +p_xl-2 # convexity
                    +2* p_xl # the lower and upper bounds for the y which we need for the new formulation
                    ,p_vl)        #Creates a new linear program model object.Specifies the number of
    #constraints (nrow) and the number of decision variables (ncol).

    #mylp<-make.lp(4*zpl,vl)
    #row.add.mode(mylp,state="on")

    set.objfn(p_mylp,obj=rep(1,p_fl),indices=(p_xl+1):(p_xl+p_fl))           #Sets the objective function. Specifies: the linear program model object,
    #the coefficients of the objective function (in this case a vector with as many 1s as flags
    #and 0s elsewhere).

    lp.control(p_mylp,sense="max")                                           #Specifies that it is a maximization problem.

    set.type(p_mylp,(p_xl+1):(p_xl+p_fl),type="binary")                      #Sets the type of the decision variables to "binary" (= "integer" with an upper
    #bound of 1 and lower bound of 0). In our case the decision variable can assume the value of 1 or 0.

    #Sets the type of the constraints, which are all less or equal ("<=")
    set.constr.type(p_mylp, c(rep(1,2*p_xl+p_xl-1+p_xl-2+2* p_xl)))

    p_iecount<-1
    p_rhs<-NULL

    #Here we compile the matrix of our linear program model object. We consider x binary, L lower bound (lower), U upper bound (upper).
    # additional definitions to not make too many mistakes with the indices
    i_x     <- 0
    i_fl    <- p_xl
    i_zl    <- p_xl+p_fl

    #P_zpl:= f_i*x_(i+1). Flag i times put price i+1
    # z<=y-L(1-f)
    # z>=y-U(1-f)

    # z-y-Lf<=-L
    # -z+Uf+y<=U


    # f1 x2
    # first two constraints with flagbar, second two with flags
    for (i in 1:p_xl){

      set.row(p_mylp,p_iecount, xt=c(1,-1,G),indices=c(i_zl+i, i_x+i, i_fl+i))    #Third constraint corresponding to y+Ux-z
      p_iecount <- p_iecount + 1

      set.row(p_mylp, p_iecount,xt=c(-1,1,G), indices=c(i_zl+i,i_x+i,i_fl+i))  #Fourth constraint corresponding to z-Lx-y
      p_iecount <- p_iecount + 1

      p_rhs<-c(p_rhs,c(G,G))                                     #The right side of the constraints corresponding to the vector [0,0, U, -L]
    }

    # are options increasing in x? x_i - x_(i+1) <= 0.
    # f1 (x2-x1)

    for (i in 2:(p_xl)){
      set.row(p_mylp, p_iecount,xt=c(1,-1), indices=c(i_zl+i,i_zl+i-1))
      p_iecount <- p_iecount + 1

      p_rhs<-c(p_rhs,-epsim)

    }

    # are options convex in x?
    # O(K2)<=alpha O(K1)+(1-alpha) O(K3)
    # we do it with flags here
    # f2*(O(K2)-alpha O(k1)-(1-alpha)O(K3))<=0

    for (i in 2:(p_xl-1)){
      alpha<-(puts$x[i]-puts$x[i+1])/(puts$x[i-1]-puts$x[i+1])
      set.row(p_mylp, p_iecount,xt=c(1,-alpha,alpha-1),
              indices=c(i_zl+i,i_zl+i-1,i_zl+i+1))
      p_iecount <- p_iecount + 1

      p_rhs<-c(p_rhs,-epsic)

    }



    #lower and upper bounds for x_i
    for (i in 1:p_xl){

      set.row(p_mylp,p_iecount, xt=c(-1),indices=c( i_x+i))    #Third constraint corresponding to y+Ux-z
      p_iecount <- p_iecount + 1

      set.row(p_mylp, p_iecount,xt=c(1), indices=c( i_x+i))  #Fourth constraint corresponding to z-Lx-y
      p_iecount <- p_iecount + 1

      p_rhs<-c(p_rhs,c(-puts$lower[i],puts$upper[i]))                                     #The right side of the constraints corresponding to the vector [0,0, U, -L]
    }


    set.rhs(p_mylp,p_rhs)                          #sets the right hand side of the linear program model object


    # lp.control(p_mylp,bb.depthlimit=0)


    p_sol<-solve(p_mylp)


    p_flags<-as.logical(get.variables(p_mylp)[(p_xl+1):(p_xl+p_fl)])
    # print(p_flags)
    solvec_dec<-get.variables(p_mylp)[(p_xl+p_fl+1):(p_xl+p_fl+p_xl)]





    # solvec_dec = solvec_dec[-which(p_flags==FALSE)]
    if(visual){
      matplot(puts$x, cbind(puts$lower,solvec_dec,puts$upper),type='l')
    }

    return(list(p_flags, solvec_dec[p_flags],p_sol))

  }






  if (increasing == TRUE){
    return(DIconvex_I(x, lower, upper, epsim, epsic,visual))
  } else {
    return(DIconvex_D(x, lower, upper,epsim, epsic,visual))
  }


}
