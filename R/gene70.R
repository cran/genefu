`gene70` <-
function(data, logged2=TRUE, verbose=FALSE) {

	if(!logged2) { data <- log2(data) }

	gt <- nrow(sig.gene70)
	data <- data[ ,intersect(dimnames(sig.gene70)[[1]], dimnames(data)[[2]])]
	genes70 <- sig.gene70[dimnames(data)[[2]], ]
	gm <- nrow(genes70)

	if(verbose && gm != gt) { warning(sprintf("%i/%i probes are used to compute the score", gm, gt)) }

	score <- apply(X=data, MARGIN=1, FUN=cor, y=sig.gene70[ ,"average.good.prognosis.profile"], method="spearman", use="complete.obs")
	score <- -score
	official.cutoff <- -0.3
	## cutoff leaving 59% of patients in the poor prognosis group in the original dataset
	risk <- ifelse(score >= official.cutoff, 1, 0)

	names(score) <- names(risk) <- dimnames(data)[[1]]
	
	return(list("score"=score, "risk"=risk))
}
