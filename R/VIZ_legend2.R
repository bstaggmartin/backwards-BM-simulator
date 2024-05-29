# legend2<-legend
# formals(legend2)<-c(formals(legend2),fill.lwd=list(NULL),fill.lty=list(NULL))
# tmp<-body(legend2)[[42]][[3]][[2]][[3]][[3]] #don't know how to find this automatically...
# tmp$lty<-quote(fill.lty)
# tmp$lwd<-quote(fill.lwd)
# body(legend2)[[42]][[3]][[2]][[3]][[3]]<-tmp
#not sure if it's appropriate since the namespace remains graphics...

#' @export
legend2<-function (x, y = NULL, legend, fill = NULL, col = par("col"),
                   border = "black", lty, lwd, pch, angle = 45, density = NULL, 
                   bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"), 
                   box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd, 
                   xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5),
                   text.width = NULL, text.col = par("col"), text.font = NULL, 
                   merge = do.lines && has.pch, trace = FALSE, plot = TRUE, 
                   ncol = 1, horiz = FALSE, title = NULL, inset = 0, xpd, title.col = text.col[1], 
                   title.adj = 0.5, title.cex = cex[1], title.font = text.font[1], 
                   seg.len = 2,
                   fill.lwd = NULL, fill.lty = NULL, ribbon.first = TRUE) 
{
  if (missing(legend) && !missing(y) && (is.character(y) || 
                                         is.expression(y))) {
    legend <- y
    y <- NULL
  }
  mfill <- !missing(fill) || !missing(density)
  if (!missing(xpd)) {
    op <- par("xpd")
    on.exit(par(xpd = op))
    par(xpd = xpd)
  }
  text.font <- if (is.null(text.font)) 
    par("font")
  else text.font
  title <- as.graphicsAnnot(title)
  if (length(title) > 1) 
    stop("invalid 'title'")
  legend <- as.graphicsAnnot(legend)
  if (any(sapply(legend, is.language))) 
    legend <- as.expression(legend)
  n.leg <- length(legend)
  if (n.leg == 0) 
    stop("'legend' is of length 0")
  auto <- if (is.character(x)) 
    match.arg(x, c("bottomright", "bottom", "bottomleft", 
                   "left", "topleft", "top", "topright", "right", "center"))
  else NA
  if (is.na(auto)) {
    xy <- xy.coords(x, y, setLab = FALSE)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2) 
      stop("invalid coordinate lengths")
  }
  else nx <- 0
  reverse.xaxis <- par("xaxp")[1] > par("xaxp")[2]
  reverse.yaxis <- par("yaxp")[1] > par("yaxp")[2]
  xlog <- par("xlog")
  ylog <- par("ylog")
  cex <- rep(cex, length.out = n.leg)
  x.intersp <- rep(x.intersp, length.out = n.leg)
  seg.len <- rep(seg.len, length.out = n.leg)
  rect2 <- function(left, top, dx, dy, density = NULL, angle, 
                    ...) {
    r <- left + dx
    if (xlog) {
      left <- 10^left
      r <- 10^r
    }
    b <- top - dy
    if (ylog) {
      top <- 10^top
      b <- 10^b
    }
    rect(left, top, r, b, angle = angle, density = density, 
         ...)
  }
  segments2 <- function(x1, y1, dx, dy, ...) {
    x2 <- x1 + dx
    if (xlog) {
      x1 <- 10^x1
      x2 <- 10^x2
    }
    y2 <- y1 + dy
    if (ylog) {
      y1 <- 10^y1
      y2 <- 10^y2
    }
    segments(x1, y1, x2, y2, ...)
  }
  
  
  
  matplot2 <- function(x1, y1, dx, dy, ...){
    nn <- max(4, round(dx))
    a <- rep(c(-1, 0), length.out = nn)
    b <- rep(c(0, 1), length.out = nn)
    xx <- do.call(cbind,lapply(seq_along(x1), function(ii) seq(x1[ii], x1[ii] + dx[ii], length.out = nn)))
    yy <- matrix(y1, nn, length(y1), byrow = TRUE) + 0.5 * yc * runif(nn * length(y1), a, b)
    matplot(x = xx, y = yy, ..., add = TRUE, type = "l")
  }
  
  
  
  points2 <- function(x, y, ...) {
    if (xlog) 
      x <- 10^x
    if (ylog) 
      y <- 10^y
    points(x, y, ...)
  }
  text2 <- function(x, y, ...) {
    if (xlog) 
      x <- 10^x
    if (ylog) 
      y <- 10^y
    text(x, y, ...)
  }
  colwise <- function(x, n, ncol, n.legpercol, fun, reverse = FALSE) {
    xmat <- matrix(c(rep(x, length.out = n), rep(0L, n.legpercol * 
                                                   ncol - n)), ncol = ncol)
    res <- apply(xmat, 2, fun)
    res[res == 0L] <- max(res)
    if (reverse) 
      -res
    else res
  }
  rowwise <- function(x, n, ncol, n.legpercol, fun, reverse = FALSE) {
    xmat <- matrix(c(rep(x, length.out = n), rep(0L, n.legpercol * 
                                                   ncol - n)), ncol = ncol)
    res <- apply(xmat, 1, fun)
    if (reverse) 
      -res
    else res
  }
  if (trace) {
    catn <- function(...) do.call(cat, c(lapply(list(...), 
                                                formatC), "\n"))
    fv <- function(...) paste(vapply(lapply(list(...), formatC), 
                                     paste, collapse = ",", ""), collapse = ", ")
  }
  n.legpercol <- if (horiz) {
    if (ncol != 1) 
      warning(gettextf("horizontal specification overrides: Number of columns := %d", 
                       n.leg), domain = NA)
    ncol <- n.leg
    1
  }
  else ceiling(n.leg/ncol)
  Cex <- cex * par("cex")
  if (is.null(text.width)) 
    text.width <- max(abs(mapply(strwidth, legend, cex = cex, 
                                 font = text.font, MoreArgs = list(units = "user"))))
  else if ((length(text.width) > 1L && any(is.na(text.width))) || 
           (all(!is.na(text.width)) && (!is.numeric(text.width) || 
                                        any(text.width < 0)))) 
    stop("'text.width' must be numeric, >= 0, or a scalar NA")
  if (auto.text.width <- all(is.na(text.width))) {
    text.width <- abs(mapply(strwidth, legend, cex = cex, 
                             font = text.font, MoreArgs = list(units = "user")))
    ncol <- ceiling(n.leg/n.legpercol)
  }
  xyc <- xyinch(par("cin"), warn.log = FALSE)
  xc <- Cex * xyc[1L]
  yc <- Cex * xyc[2L]
  if (any(xc < 0)) 
    text.width <- -text.width
  xchar <- xc
  xextra <- 0
  y.intersp <- rep(y.intersp, length.out = n.legpercol)
  yextra <- rowwise(yc, n = n.leg, ncol = ncol, n.legpercol = n.legpercol, 
                    fun = function(x) max(abs(x)), reverse = reverse.yaxis) * 
    (y.intersp - 1)
  ymax <- sign(yc[1]) * max(abs(yc)) * max(1, mapply(strheight, 
                                                     legend, cex = cex, font = text.font, MoreArgs = list(units = "user"))/yc)
  ychar <- yextra + ymax
  ymaxtitle <- title.cex * par("cex") * xyc[2L] * max(1, strheight(title, 
                                                                   cex = title.cex, font = title.font, units = "user")/(title.cex * 
                                                                                                                          par("cex") * xyc[2L]))
  ychartitle <- yextra[1] + ymaxtitle
  if (trace) 
    catn("  xchar=", fv(xchar), "; (yextra, ychar)=", fv(yextra, 
                                                         ychar))
  if (mfill) {
    xbox <- xc * 0.8
    ybox <- yc * 0.5
    dx.fill <- max(xbox)
  }
  do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 
                                                            0))) || !missing(lwd)
  has.pch <- !missing(pch) && length(pch) > 0
  if (do.lines) {
    x.off <- if (merge) 
      -0.7
    else 0
  }
  else if (merge) 
    warning("'merge = TRUE' has no effect when no line segments are drawn")
  if (has.pch) {
    if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L], 
                                                      type = "c") > 1) {
      if (length(pch) > 1) 
        warning("not using pch[2..] since pch[1L] has multiple chars")
      np <- nchar(pch[1L], type = "c")
      pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
    }
    if (!is.character(pch)) 
      pch <- as.integer(pch)
  }
  if (is.na(auto)) {
    if (xlog) 
      x <- log10(x)
    if (ylog) 
      y <- log10(y)
  }
  if (nx == 2) {
    x <- sort(x)
    y <- sort(y)
    left <- x[1L]
    top <- y[2L]
    w <- diff(x)
    h <- diff(y)
    w0 <- w/ncol
    x <- mean(x)
    y <- mean(y)
    if (missing(xjust)) 
      xjust <- 0.5
    if (missing(yjust)) 
      yjust <- 0.5
  }
  else {
    yc <- rowwise(yc, n.leg, ncol, n.legpercol, fun = function(x) max(abs(x)), 
                  reverse = reverse.yaxis)
    h <- sum(ychar) + yc[length(yc)] + (!is.null(title)) * 
      ychartitle
    xch1 <- colwise(xchar, n.leg, ncol, n.legpercol, fun = function(x) max(abs(x)), 
                    reverse = reverse.xaxis)
    x.interspCol <- colwise(x.intersp, n.leg, ncol, n.legpercol, 
                            fun = max)
    seg.lenCol <- colwise(seg.len, n.leg, ncol, n.legpercol, 
                          fun = max)
    text.width <- colwise(text.width, n = if (auto.text.width) 
      n.leg
      else ncol, ncol, n.legpercol = if (auto.text.width) 
        n.legpercol
      else 1, fun = function(x) max(abs(x)), reverse = reverse.xaxis)
    w0 <- text.width + (x.interspCol + 1) * xch1
    if (mfill) 
      w0 <- w0 + dx.fill
    if (do.lines) 
      w0 <- w0 + (seg.lenCol + x.off) * xch1
    w <- sum(w0) + 0.5 * xch1[ncol]
    if (!is.null(title) && (abs(tw <- strwidth(title, units = "user", 
                                               cex = title.cex, font = title.font) + 0.5 * title.cex * 
                                par("cex") * xyc[1L])) > abs(w)) {
      xextra <- (tw - w)/2
      w <- tw
    }
    if (is.na(auto)) {
      left <- x - xjust * w
      top <- y + (1 - yjust) * h
    }
    else {
      usr <- par("usr")
      inset <- rep_len(inset, 2)
      insetx <- inset[1L] * (usr[2L] - usr[1L])
      left <- switch(auto, bottomright = , topright = , 
                     right = usr[2L] - w - insetx, bottomleft = , 
                     left = , topleft = usr[1L] + insetx, bottom = , 
                     top = , center = (usr[1L] + usr[2L] - w)/2)
      insety <- inset[2L] * (usr[4L] - usr[3L])
      top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] + 
                      h + insety, topleft = , top = , topright = usr[4L] - 
                      insety, left = , right = , center = (usr[3L] + 
                                                             usr[4L] + h)/2)
    }
  }
  if (plot && bty != "n") {
    if (trace) 
      catn("  rect2(", left, ",", top, ", w=", w, ", h=", 
           h, ", ...)", sep = "")
    rect2(left, top, dx = w, dy = h, col = bg, density = NULL, 
          lwd = box.lwd, lty = box.lty, border = box.col)
  }
  xt <- left + xc + xextra + rep(c(0, cumsum(w0))[1L:ncol], 
                                 each = n.legpercol, length.out = n.leg)
  topspace <- 0.5 * ymax + (!is.null(title)) * ychartitle
  yt <- top - topspace - cumsum((c(0, ychar)/2 + c(ychar, 0)/2)[1L:n.legpercol])
  yt <- rep(yt, length.out = n.leg)
  
  
  
  if(ribbon.first){
    if (mfill) {
      if (plot) {
        if (!is.null(fill)) 
          fill <- rep_len(fill, n.leg)
        rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox, 
              col = fill, density = density, angle = angle, 
              border = border,
              lwd = fill.lwd,
              lty = fill.lty)
      }
      xt <- xt + dx.fill
    }
  }else if(mfill){
    if (plot){
      xt <- xt + dx.fill
    }
  }
  
  
  
  if (plot && (has.pch || do.lines)) 
    col <- rep_len(col, n.leg)
  if (missing(lwd) || is.null(lwd)) 
    lwd <- par("lwd")
  
  
  
  if (do.lines) {
    if (missing(lty) || is.null(lty)) 
      lty <- 1
    lty <- rep_len(lty, n.leg)
    lwd <- rep_len(lwd, n.leg)
    ok.l <- !is.na(lty) & (is.character(lty) | lty > 0) & 
      !is.na(lwd)
    if (trace) 
      catn("  matplot2(", xt[ok.l] + x.off * xchar, ",",
           yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
      # catn("  segments2(", xt[ok.l] + x.off * xchar, ",", 
      #      yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
    if (plot) 
      matplot2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len * xchar,
               dy = 0, lty = lty[ok.l], lwd = lwd[ok.l],
               col = col[ok.l])
      # segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len * 
      #             xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l], 
      #           col = col[ok.l])
    xt <- xt + (seg.len + x.off) * xchar
  }
  if(!ribbon.first){
    if (mfill) {
      if (plot) {
        if (!is.null(fill)) 
          fill <- rep_len(fill, n.leg)
        tmp.xt <- if(merge) xt + x.off * xchar else xt - (x.off + seg.len) * xchar
        rect2(left = tmp.xt, top = yt + ybox/2, dx = xbox, dy = ybox, 
              col = fill, density = density, angle = angle, 
              border = border,
              lwd = fill.lwd,
              lty = fill.lty)
      }
    }
  }
  
  
  
  if (has.pch) {
    pch <- rep_len(pch, n.leg)
    pt.bg <- rep_len(pt.bg, n.leg)
    pt.cex <- rep_len(pt.cex, n.leg)
    pt.lwd <- rep_len(pt.lwd, n.leg)
    ok <- !is.na(pch)
    if (!is.character(pch)) {
      ok <- ok & (pch >= 0 | pch <= -32)
    }
    else {
      ok <- ok & nzchar(pch)
    }
    x1 <- (if (merge && do.lines) 
      xt - (seg.len/2) * xchar
      else xt)[ok]
    y1 <- yt[ok]
    if (trace) 
      catn("  points2(", x1, ",", y1, ", pch=", pch[ok], 
           ", ...)")
    if (plot) 
      points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok], 
              bg = pt.bg[ok], lwd = pt.lwd[ok])
  }
  xt <- xt + x.intersp * xc
  if (plot) {
    if (!is.null(title)) 
      text2(left + w * title.adj, top - ymaxtitle, labels = title, 
            adj = c(title.adj, 0), cex = title.cex, col = title.col, 
            font = title.font)
    text2(xt, yt, labels = legend, adj = adj, cex = cex, 
          col = text.col, font = text.font)
  }
  invisible(list(rect = list(w = w, h = h, left = left, top = top), 
                 text = list(x = xt, y = yt)))
}
