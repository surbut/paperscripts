effect_plot=
function (model, pred, pred.values = NULL, centered = "all", 
          plot.points = FALSE, interval = FALSE, data = NULL, at = NULL, 
          int.type = c("confidence", "prediction"), int.width = 0.95, 
          outcome.scale = "response", robust = FALSE, cluster = NULL, 
          vcov = NULL, set.offset = 1, x.label = NULL, y.label = NULL, 
          pred.labels = NULL, main.title = NULL, colors = "black", 
          line.colors = colors, line.thickness = 1.1, point.size = 1.5, 
          point.alpha = 0.6, jitter = 0, rug = FALSE, rug.sides = "lb", 
          force.cat = FALSE, cat.geom = c("point", "line", "bar"), 
          cat.interval.geom = c("errorbar", "linerange"), cat.pred.point.size = 3.5, 
          partial.residuals = FALSE, color.class = colors, ...) 
{
  pred <- quo_name(enexpr(pred))
  if ("interval" %nin% names(match.call())[-1] && !(is.numeric(get_data(model, 
                                                                        warn = FALSE)[[pred]]) && force.cat == FALSE)) {
    interval <- TRUE
  }
  if (force.cat == TRUE && is.null(pred.values)) {
    if (is.null(data)) {
      data <- get_data(model)
    }
    pred.values <- sort(unique(suppressMessages(data[[pred]])))
  }
  if (!all(color.class == colors)) 
    colors <- color.class
  colors <- get_colors(colors)
  pred_out <- make_predictions(model, pred = pred, pred.values = pred.values, 
                               at = at, center = centered, interval = interval, int.type = int.type, 
                               int.width = int.width, outcome.scale = outcome.scale, 
                               robust = robust, cluster = cluster, vcov = vcov, set.offset = set.offset, 
                               return.orig.data = TRUE, partial.residuals = partial.residuals, 
                               data = data, ...)
  pm <- pred_out[[1]]
  d <- pred_out[[2]]
  dots <- list(...)
  if (!is.null(dots$dpar) && plot.points == TRUE) {
    plot.points <- FALSE
    warn_wrap("The plot.points argument is not compatible with distributional\n              parameters specified in `dpar`.")
  }
  if (!is.null(dots$point.color)) {
    warn_wrap("The 'point.color' argument is deprecated and is now ignored. \n               You can change the color of points with the 'colors' argument.")
  }
  if (is.numeric(d[[pred]]) && force.cat == FALSE) {
    plot_effect_continuous(predictions = pm, pred = pred, 
                           plot.points = plot.points | partial.residuals, interval = interval, 
                           data = d, x.label = x.label, y.label = y.label, pred.labels = pred.labels, 
                           main.title = main.title, colors = line.colors, line.thickness = line.thickness, 
                           jitter = jitter, resp = get_response_name(model, 
                                                                     ...), weights = get_weights(model, d)$weights_name, 
                           rug = rug, rug.sides = rug.sides, point.size = point.size, 
                           point.alpha = point.alpha, point.color = colors)
  }
  else {
    plot_cat(predictions = pm, pred = pred, data = d, geom = cat.geom, 
             pred.values = pred.values, interval = interval, plot.points = plot.points | 
               partial.residuals, pred.labels = pred.labels, 
             x.label = x.label, y.label = y.label, main.title = main.title, 
             colors = line.colors, weights = get_weights(model, 
                                                         d)$weights_name, resp = get_response_name(model, 
                                                                                                   ...), jitter = jitter, interval.geom = cat.interval.geom, 
             line.thickness = line.thickness, point.size = point.size, 
             pred.point.size = cat.pred.point.size, point.alpha = point.alpha, 
             point.color = colors)
  }
}