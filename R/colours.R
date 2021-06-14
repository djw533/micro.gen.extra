random_x_colours <- function(num,hex) {
  color = grDevices::colors()[grep('gr(a|e)y|white', grDevices::colors(), invert = T)]
  col_vector <- sample(color, num)

  #set the hex color strings (if specified) - otherwise return strings:
  if (isTRUE(hex)) {
    hex_colours <- list()
    for (col in col_vector) {
      hex_colours <- append(hex_colours,gplots::col2hex(col))
      hex_colours <- as.character(hex_colours)
    }
    return(hex_colours)
  } else {
    return(col_vector)
  }
}


color.bar <- function(lut, min, max, axis_bool, alpha_num, ylab_string, title_string, nticks=8, ticks=seq(min, max, len=nticks)) {
  scale = (length(lut)-1)/(max-min)
  #dev.new(width=1.75, height=5)
  if (axis_bool == T) {
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', yaxt='n', xlab='', main = title_string, ylab=ylab_string)
    axis(2, ticks, las=1)
  } else {
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', yaxt='n', xlab='', ylab='', main = title_string)
  }
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=alpha(lut[i], alpha_num), border=NA)
  }
}
