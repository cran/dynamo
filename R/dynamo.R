.packageName <- "dynamo"

.First.lib <- function(lib, pkg)
{
     library.dynam("dynamo", pkg, lib)
}

`dm.spec` <- 
function( frml, innovations, nobs=NULL, data)
{
	# utilities
	getexpl <- function( args, frm)
	{
		# get variables names
		vars <- attr( terms( formula(paste('~',args))) , 'term.labels' )
		# get variables from frame and put them in a matrix
		expl <- numeric()
		for( var in vars )
		{
			expl <- cbind( expl , frm[[var]] )
		}

		expl
	}

	parse.args <- function( args.type, args, frm)
	{
		res <- c()

		for( type in args.type )
		{
			if( type == 'i' )
			{
				if( regexpr('[0-9]',args[1])[[1]]==-1 ) break;
			}

			if( type == 'x' )
			{
			}
		}
	}

	# map
	spec.map.mean <- list( none='DYNAMO_MEAN_NONE', const='DYNAMO_MEAN_CONST', arma='DYNAMO_MEAN_ARMA', carma='DYNAMO_MEAN_CARMA', mem='DYNAMO_MEAN_MEM')
	spec.map.var <- list( none='DYNAMO_VAR_NONE', const='DYNAMO_VAR_CONST', garch='DYNAMO_VAR_GARCH')
	spec.map.news <- list( norm='DYNAMO_NEWS_NORM', st='DYNAMO_NEWS_ST', sskt='DYNAMO_NEWS_SSKT', exp='DYNAMO_NEWS_EXP', gamma='DYNAMO_NEWS_GAMMA', weibull='DYNAMO_NEWS_WEIBULL')

	# labels
	labels.mean = c('ar','ma','arma','mem','acd','smem')
	labels.var = c('arch','garch','gjrgarch','egarch','aparch','sgarch')
	labels.all = c(labels.mean,labels.var)
 
	# variable initialization
	frml.char <- as.character(frml) 

	if( length(frml.char)==2 )
	{ 
		y <- rep(0,nobs);
	}
	else
	{ 
		y <- data[[ frml.char[2] ]]
		nobs = length(y)
	}

	X = matrix( 0 , 0 , 0 )
	Z = matrix( 0 , 0 , 0 )

	# IMPORTANT: substitute all interaction term '*' between dynamo terms with '%*%'
	if( length(frml.char)==3 ) frml <- formula( paste( frml.char[2] , '~' , sub(') \\*',')%*%',frml.char[3])) )
	else frml <- formula( paste( '~' , sub(') \\*',')%*%',frml.char[2])) )

	# get terms
	intercept  <- attr( terms(frml), 'intercept')
	frml.terms <- attr( terms(frml), 'term.labels')

	# get terms names
	term.names <- c()
	for( str in as.character(frml.terms) )
	{
		if( regexpr("\\(",str)[[1]] != -1 )
		{
			str <- gsub( '\\(.*\\) %\\*% ' , '', str)
			str <- gsub( '\\(.*\\)' , '', str)
			term.names <- c( term.names , str)
		}
		else stop( paste("Invalid Formula Specification: all dynamo terms must have arguments", str), call.=FALSE)
	}
	term.names <- tolower( term.names )

	# get terms args
	term.args <- c()
	for( i in 1:length(frml.terms) )
	{
		str <- frml.terms[i]

		str <- gsub( '\\) %\\*% .*\\(' , ',', str)
		str <- gsub( '.*\\(' , '', str)
		str <- gsub( '\\).*' , '', str)

		term.args[[i]] <- strsplit(str,',')[[1]]
	}

	# some boring checks...
	# check if all terms are known
	if( sum( (term.names %in% labels.all)==FALSE )>0 )
	{
		stop( paste( 
			"Invalid Formula Specification: unknown term(s)",
			paste(term.names[(term.names %in% labels.all)==FALSE],collapse=', ') ))
	}

	# check if there are not term duplicates

	# check if there are not too many mean terms
	if( sum( labels.mean %in% term.names) > 1 ) 
	{
		stop(paste(
		"Invalid Formula Specification: too many mean terms",
		paste( labels.mean[(labels.mean %in% term.names)] , collapse=', ' ) ))
	}
	# check if there are not too many var terms
	if( sum( labels.var %in% term.names) > 1 )
	{
		stop(paste(
		"Invalid Formula Specification: too many variance terms",
		paste( labels.var[(labels.var %in% term.names)] , collapse=', ' ) ))
	}
	# check if terms are in the right order
	if( length(term.names)>1 )
	{
		m <- diff( match( term.names , labels.all ) )
		if( min(m)<0 )
		{
			stop('Invalid Formula Specification: bad term ordering')
		}
	}

	# specification reconstruction loop 
	attr = c(spec.map.mean[['none']],spec.map.var[['none']],spec.map.news[['norm']])
	num.mean = c(0,0,0,0) 
	num.var = c(0,0,0,0) 
	valid=FALSE
	
	# set inovations distribution (if valid)
	if( !is.null(innovations) )
	{
		if( !(innovations %in% names(spec.map.news)) ) stop('Invalid Formula Specification: bad innovation distribution')	

		attr[3] = spec.map.news[[innovations]]
	}
	# default innovation distribution (depending on class of models)
	else
	{
		if( term.names[1] %in% c('mem','acd','smem') )
		{
			attr[3] = spec.map.news[['exp']]
		}
		else 
		{
			attr[3] = spec.map.news[['norm']]
		}
	}

	for( i in 1:length(frml.terms) )
	{
		#print( term.names[i] )

		if( term.names[i] == 'ar' )
		{
			if( intercept ) attr[1] = spec.map.mean[['carma']]
			else attr[1] = spec.map.mean[['arma']]

			valid=TRUE
		}
		if( term.names[i] == 'ma' )
		{
			if( intercept ) attr[1] = spec.map.mean[['carma']]
			else attr[1] = spec.map.mean[['arma']]

			valid=TRUE
		}
		if( term.names[i] == 'arma' )
		{
			if( intercept ) attr[1] = spec.map.mean[['carma']]
			else attr[1] = spec.map.mean[['arma']]

			args = term.args[[i]]

			# checks
			#if( length(args)!=2 && length(args)!=3 ) break;
			#if( regexpr('[0-9]',args[1])[[1]]==-1 || regexpr('[0-9]',args[2])[[1]]==-1 ) break;

			#num.mean[1] = as.integer( args[1] )
			#num.mean[2] = as.integer( args[2] )

			#if( length(args)==3 )
			#{
			#	X = getexpl(args[3],parent.frame());
			#	num.mean[3] = ncol(X)
			#}

			valid=TRUE
		}
		if( term.names[i] == 'mem' || term.names[i] == 'acd' )
		{
			attr[1] = spec.map.mean[['mem']]

			args = term.args[[i]]

			# checks
			#if( length(args)!=2 && length(args)!=3 ) break;
			#if( regexpr('[0-9]',args[1])[[1]]==-1 || regexpr('[0-9]',args[2])[[1]]==-1 ) break;

			num.mean[1] = as.integer( args[1] )
			num.mean[2] = as.integer( args[2] )

			if( length(args)==3 ) 
			{
				X = getexpl( args[3], data)
				num.mean[3] = ncol(X)
			}

			#ret <- parse.args( c('i','i','x') , args, data);
			#num.mean <- ret$num;
			#X <- ret

			valid=TRUE
		}
		if( term.names[i] == 'smem' || term.names[i] == 'sacd' )
		{
			valid=TRUE
		}
		if( term.names[i] == 'arch' )
		{
			attr[2] = spec.map.var[['garch']]

			args = term.args[[i]]

			num.var[1] = as.integer( 0 )
			num.var[2] = as.integer( args[1] )

			if( length(args)==2 ) 
			{
				Z = getexpl( args[2], data)
				num.var[3] = ncol(Z)
			}

			valid=TRUE
		}
		if( term.names[i] == 'garch' )
		{
			attr[2] = spec.map.var[['garch']]

			args = term.args[[i]]

			num.var[1] = as.integer( args[1] )
			num.var[2] = as.integer( args[2] )

			if( length(args)==3 ) 
			{
				Z = getexpl( args[3], data)
				num.var[3] = ncol(Z)
			}

			valid=TRUE
		}
		if( term.names[i] == 'gjrgarch' )
		{
			attr[2] = spec.map.var[['gjrgarch']]

			#num.var[1] = as.integer( args[1] )
			#num.var[2] = as.integer( args[2] )

			valid=TRUE
		}
	}

	if( !valid ) stop("Invalid Formula Specification")

	list( attr=attr, num.mean=num.mean, num.var=num.var, nobs=nobs, y=y , X=X , Z=Z )
}

`dm` <-
function(formula,innovations=NULL,data=parent.frame(),est='mle',maxiter=NULL,param=NULL,ftol=NULL,log='none')
{
	# maps
	opt.map.obj <- list( mle="DYNAMO_OBJ_LOGLIK", ridge="DYNAMO_OBJ_PLOGLIK_RIDGE", gridge="DYNAMO_OBJ_PLOGLIK_GRIDGE");
	opt.map.loglev <- list( none="DYNAMO_LOGLEV_NONE", basic="DYNAMO_LOGLEV_BASIC", detail="DYNAMO_LOGLEV_DETAIL");
	
	cl <- match.call()

	# get specification
	spec <- dm.spec( formula, innovations, NULL, data)	

	# get options
	if( is.null(maxiter) ){ maxiter <- 0 }
	if( is.null(ftol) ){ ftol <- 0 }
	if( is.null(param) ){ init.param <- 0 } else { init.param <- 1 }
	opt <- list( obj=opt.map.obj[est], loglev=opt.map.loglev[log], maxiter=maxiter, ftol=ftol, init.param=init.param );

	if( is.null(param) ) param <- rep(0,1000)

	# call routine
	res <- .C("dynamo_fit_wrapper",
		info = as.integer(0),
		info.msg = as.character( paste( rep(" ",1000) , collapse='') ),
		y = as.vector( spec$y ),
		cond = as.vector(rep(0,4*spec$nobs)), # there can be no more than *4* cond components
		resid = as.vector( spec$y ),
		nobs = as.integer(spec$nobs),
		ncnd = as.integer(0),
		spec.attr = as.character(spec$attr),
		spec.num = as.integer( c(spec$num.mean,spec$num.var) ),
		param = as.vector(param),
		avc = as.vector( spec$y ),
		gradient = as.vector(param),
		obj = as.vector(1),
		nparam = as.integer(length(param)),
		spec$X,
		spec$Z,
		fitattr = as.character( c(opt$obj,opt$loglev) ),
		fitnum = as.integer( c(opt$maxiter,opt$init.param,opt$ftol) ),
		PACKAGE = 'dynamo' )

	if( res$info ){ stop( res$info.msg, call.=FALSE) }

	obj <- list( call=cl, spec=spec, coefficients=res$param[1:res$nparam], vc=matrix(res$avc[1:(res$nparam**2)],res$nparam,res$nparam), residuals=res$resid, gradient=res$gradient[1:res$nparam], cond=t(matrix(res$cond,res$ncnd,res$nobs)), obj=res$obj, y=spec$y, nobs=res$nobs, ncnd=res$ncnd )
	obj$mu <- obj$cond[,1]
	obj$sigma2 <- obj$cond[,2]
	obj$aic <- -2*obj$obj-2*length(obj$coefficients);
	obj$bic <- -2*obj$obj-length(obj$coefficients)*log(length(obj$y))

	class(obj) <- 'dm'

	obj
}

`rdm` <-
function (formula,innovations=NULL,n,param,data=parent.frame(),seed=as.integer(Sys.time()),lst=FALSE)
{
	cl <- match.call()

	spec <- dm.spec( formula, innovations, n, data)

	res <- .C("dynamo_sim_wrapper",
		info = as.integer(0),
		info.msg = as.character( paste( rep(" ",1000) , collapse='') ),
		y = as.vector( spec$y ),
		cond = as.vector(rep(0,4*n)), # there can be no more than *4* cond components
		nobs = as.integer(n),
		ncnd = as.integer(0),
		spec.attr = as.character(spec$attr),
		spec.num = as.integer( c(spec$num.mean,spec$num.var) ),
		param = as.vector(param),
		nparam = as.integer(length(param)),
		spec$X,
		spec$Z,
		seed = as.integer(seed), 
		PACKAGE = 'dynamo' )

	if( res$info ){ stop( res$info.msg, call.=FALSE) }

	if( !lst )
	{
		out <- res$y
	}
	else
	{
		cond <- t( matrix( res$cond, res$ncnd, res$nobs))
		out <- list( y=res$y, mu=cond[,1], sig2=cond[,2], call=cl)
	}

	out
}

print.dm <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

plot.dm <- function( x , ... )
{
	plot( x$y , col='dark grey' )
	lines( x$mu , col='blue' )
}

summary.dm <- function( object , diag.resid.lag=5 , ... )
{
	s.obj <- list( dm=object )

	s.obj$Qz  <- Box.test( object$resid   , lag = diag.resid.lag , type = c("Ljung-Box") )
	s.obj$Qz2 <- Box.test( object$resid^2 , lag = diag.resid.lag , type = c("Ljung-Box") )

	s.obj$mse <- sum( (object$y-object$mu)^2 )/length(object$y)
	s.obj$mae <- sum( abs(object$y-object$mu) )/length(object$y)

	s.obj$aic <- 2*object$obj+2*length(object$coef);
	s.obj$bic <- 2*object$obj+length(object$coef)*log(length(object$spec$y));

	class( s.obj ) <- 'summary.dm'

	s.obj
}

print.summary.dm <- function( x , ... )
{
	cat("\nCall:\n", deparse(x$dm$call), "\n\n", sep = "")
}

predict.dm <- function( object, type=stop('Forecast type not specified!',call.=FALSE), iterate=NULL, nfor=1, ql=0.05, qu=0.95, yo=NULL, Xo=matrix(0), Zo=matrix(0), ... )
{
	pred.map.obj <- list( dynamic="DYNAMO_PRED_DYNAMIC", static='DYNAMO_PRED_STATIC' );

	if( !is.null(yo) ) nfor = length(yo)+1;

	res <- .C("dynamo_pred_wrapper",
		info = as.integer(0),
		info.msg = as.character( paste( rep(" ",1000) , collapse='') ),
		pred = as.character( pred.map.obj[type]  ),
		quant = as.vector( c(ql,qu) ),
		cond_pred = as.vector(rep(0,4*nfor)), # there can be no more than *4* cond components
		quant_pred = as.vector(rep(0,2*nfor)),
		spec.attr = as.character(object$spec$attr),
		spec.num = as.integer( c(object$spec$num.mean,object$spec$num.var) ),
		param = as.vector(object$coef),
		y = as.vector( object$spec$y ),
		X=object$spec$X,
		Z=object$spec$Z,
		nobs = as.integer(object$spec$nobs),
		yo=as.vector(yo),
		Xo=as.matrix(Xo),
		Zo=as.matrix(Zo),
		nfor=as.integer(nfor), 
		PACKAGE = 'dynamo' )

	out <- list( cond_pred=t(matrix(res$cond_pred,object$ncnd,nfor)) , quant_pred=t(matrix(res$quant_pred,2,nfor)) ) 
	out$mu <- out$cond_pred[,1]
	out$sigma2 <- out$cond_pred[,2]

	out
}

