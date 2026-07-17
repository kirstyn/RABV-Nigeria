function(tree_view = NULL, node, xexpand = NULL) {
  p <- ggtree:::get_tree_view(tree_view)
  sp <- ggtree:::offspring(p, node, self_include = TRUE)
  
  xr <- range(sp$x)
  
  if (is.null(xexpand)) {
    xr[2] <- xr[2] + diff(xr)/10
  } else {
    xr[2] <- xr[2] + diff(xr)/10 + xr[2] * xexpand
  }
  
  yr <- range(sp$y)
  
  facet_zoom <- getFromNamespace("facet_zoom", "ggforce")
  
  p + facet_zoom(
    ylim = yr,
    show.area = TRUE,
    horizontal = TRUE,
  )
}
