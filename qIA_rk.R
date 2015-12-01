# RKWard plugin for the analysis of Quantitative Isothermal Amplification reactions.

require(rkwarddev)
local({
  
  # Author names and contact information
  about.info <- rk.XML.about(
    name = "qIA analysis",
    author = c(
      person(given = "Stefan", family = "Roediger",
             email = "stefan.roediger@b-tu.de", 
             role = c("aut","cre"))),
    about = list(desc = "GUI interface to analyze a quantitative isothermal amplification",
                 version = "0.0.1", url = "")
  )
  
  ## help page
  plugin.summary <- rk.rkh.summary(
    "Analysis of quantitative isothermal amplification (qIA) curve data. The plugin is primarily targeted at the analysis of qIA data from nucleic acid experiments but might be used in other cases too."
  )
  
  plugin.usage <- rk.rkh.usage(
    "Chose a data set and method for the qIA."
  )
  
  # Define dependencies
  dependencies.info <- rk.XML.dependencies(dependencies = list(rkward.min = "0.6.3"), 
                                           package = list(c(name = "chipPCR", min = "0.0.8.3")))
  ## General settings
  
  # Definition of plot labels and appearence
  generic.plot.options <- rk.plotOptions()
  plot.main <- rk.XML.input(label = "Main title", initial = "qIA")
  plot.xlab <- rk.XML.input(label = "Abscissa label", initial = "Time")
  plot.ylab <- rk.XML.input(label = "Ordinate label", initial = "Signal")
  label.frame <- rk.XML.frame(rk.XML.row(plot.main),
                              rk.XML.row(plot.xlab, plot.ylab),
                              rk.XML.stretch(), label = "Labels")
  
  var.select <- rk.XML.varselector(label = "Select data")
  
  selected.x <- rk.XML.varslot(label = "Time", source = var.select, multi = TRUE, 
                               types = "number", required = TRUE)
  
  single.x.chk <- rk.XML.cbox(label = "Single time column", value = "1", un.value = "0")
  scale.time.spin <- rk.XML.spinbox(label = "Time factor", min = "0", initial = "1", precision = 4)
  
  var.data <- rk.XML.varslot(label = "Signal", source = var.select, multi = TRUE, 
                             types = "number", required = TRUE)
  
  scale.axes.drop <- rk.XML.dropdown(label = "Scale axes",
                                     options = list("None" = c(val = "\n", chk = TRUE),
                                                    "Scale x" = c(val = "x <- log10(x)\nis.na(x) <- sapply(x, is.infinite)\n"),
                                                    "Scale y" = c(val = "y <- log10(y)\nis.na(y) <- sapply(y, is.infinite)\n"),
                                                    "Scale xy" = c(val = "x <- log10(x)\ny <- log10(y)\nis.na(x) <- sapply(x, is.infinite)\nis.na(y) <- sapply(y, is.infinite)\n")))
  
  pointtype.drop <- rk.XML.dropdown(label = "Type of all points or lines",
                                    options = list("Plot individual points" = c(val = "p"),
                                                   "Plot lines" = c(val = "l", chk = TRUE),
                                                   "Plot points connected by lines (both)" = c(val = "b"),
                                                   "Plot points overlaid by lines " = c(val = "o"),
                                                   "Plot vertical lines from points to the zero axis (high-density)" = c(val = "h"),
                                                   "Step-function plots, the left edge defines the point" = c(val = "s"),
                                                   "Step-function plots, the right edge defines the point" = c(val = "S")
                                    ))
  
  pch.spin <- rk.XML.spinbox(label = "Plotting symbol", min = "1", max = "100", initial = "20", real = FALSE)
  
  # Plot preview
  preview.chk <- rk.XML.preview(label = "Preview")
  suppress.warnings.chk <- rk.XML.cbox(label = "Show warnings", value = "0", un.value = "-1")
  generic.plot.options <- rk.plotOptions()
  plot.text.color <- rk.plotOptions(embed = "rkward::color_chooser", button = FALSE)
  
  basic.settings <- rk.XML.row(
    var.select,
    rk.XML.col(
      rk.XML.row(selected.x,
                 var.data),single.x.chk,
		 preview.chk,
		 suppress.warnings.chk,
      rk.XML.stretch()
    ))
  
  # Defintion of setting for the analysis
  # Definion for the CPP function
  smoother.chk <- rk.XML.cbox(label = "Do not use smoother", value = "1", un.value = "0")
  method.drop <- rk.XML.dropdown(label = "Smoothing method",
                                 options = list("LOWESS" = c(val = "lowess"), 
                                                "Moving average" = c(val = "mova"),
                                                "Savitzky-Golay smoothing filter" = c(val = "savgol", chk = TRUE),
                                                "Cubic spline smooth" = c(val = "smooth"),
                                                "Standard cubic spline smooth" = c(val = "spline"),
                                                "Friedman\'s SuperSmoother" = c(val = "supsmu"),
                                                "Weighted Whittaker (1st order finite difference penalty)" = c(val = "whit1"),
                                                "Weighted Whittaker (2nd order finite difference penalty)" = c(val = "whit2")))
  trans.chk <- rk.XML.cbox(label = "Correct background slope", value = "1", un.value = "0")
  method.reg.drop <- rk.XML.dropdown(label = "Method for linear regression",
                                     options = list("Linear regression" = c(val = "least"), 
                                                    "Rank-based Estimates of Regression Coefficients" = c(val = "rfit"),
                                                    "MM-type Estimators for Linear Regression" = c(val = "lmrob", chk = TRUE),
                                                    "Quantile regression" = c(val = "rq")))
  bg.outliers.chk <- rk.XML.cbox(label = "Remove outliers in the background range", value = "1", un.value = "0")
  median.chk <- rk.XML.cbox(label = "Use median instead of mean for outliers replacement", value = "1", un.value = "0")
  method.norm.drop <- rk.XML.dropdown(label = "normalization",
                                      options = list("None" = c(val = "none"), 
                                                     "Minimum-Maximum" = c(val = "minm", chk=TRUE),
                                                     "Lower and upper quantile" = c(val = "luqn"),
                                                     "Z-score" = c(val = "zscore")))
  qnL.spin <- rk.XML.spinbox(label = "qnL", min = "0.0001", max = "0.9999", initial = "0.03", precision = 3)
  amptest.chk <- rk.XML.cbox(label = "Test for a positive amplification", value = "1", un.value = "0")
  manual.chk <- rk.XML.cbox(label = "Use fixed threshold value for background", value = "1", un.value = "0")
  nL.spin <- rk.XML.spinbox(label = "Fixed threshold value for the background", min = "-1000000000", max = "1000000000", initial = "0.1", precision = 3)
  bg.range.start.spin <- rk.XML.spinbox(label = "Background start", min = "0", initial = "0", precision = 3)
  bg.range.end.spin <- rk.XML.spinbox(label = "Background end", min = "0", initial = "10", precision = 3)
  
  bg.range.frame <- rk.XML.frame(rk.XML.row(trans.chk, method.reg.drop), 
                                 rk.XML.row(bg.range.start.spin, bg.range.end.spin), 
                                 rk.XML.stretch(), label="Background correction")
  
  threshold.frame <- rk.XML.frame(rk.XML.row(manual.chk, nL.spin),
                                  rk.XML.stretch(), label="Threshold settings")
  
  smooth.frame <- rk.XML.frame(rk.XML.row(smoother.chk, method.drop),
                               rk.XML.stretch(), label="Smooth settings")
  
  normalization.frame <- rk.XML.frame(rk.XML.row(method.norm.drop, qnL.spin),
                                      rk.XML.stretch(), label="Normalization settings")
  
  
  # Quantification settings
  r.spin <- rk.XML.spinbox(label = "Fluorescence threshold", min = "-1000000000", max = "1000000000", initial = "1", precision = 3)
  auto.chk <- rk.XML.cbox(label = "Automatic estimation of the threshold", value = "1", un.value = "0")
  linear.chk <- rk.XML.cbox(label = "Quadratic regression", value = "1", un.value = "0")
  show.th.cyc.chk <- rk.XML.cbox(label = "Show Cycle threshold", value = "1", un.value = "0")
  
  quantification.frame <- rk.XML.frame(rk.XML.row(show.th.cyc.chk, r.spin),
                                       rk.XML.row(linear.chk, auto.chk),
                                       rk.XML.stretch(), label="Quantification settings")
  
  # Legend stettings
  legend.pos.drop <- rk.XML.dropdown(label = "Position of legend",
                                     options = list("Bottomright" = c(val = "bottomright"), 
                                                    "Bottom" = c(val = "bottom"),
                                                    "Bottomleft" = c(val = "bottomleft"),
                                                    "Left" = c(val = "left"),
                                                    "Topleft" = c(val = "topleft", chk = TRUE),
                                                    "Top" = c(val = "top"),
                                                    "Topright" = c(val = "topright"),
                                                    "Right" = c(val = "right"),
                                                    "Center" = c(val = "center")))
  
  ncol.legend.spin <- rk.XML.spinbox(label = "Number of columns in legend", min = "1", initial = "1", real = FALSE)
  legend.frame <- rk.XML.frame(rk.XML.row(pointtype.drop, pch.spin),
                               rk.XML.row(legend.pos.drop, ncol.legend.spin), 
                               rk.XML.stretch(), label="Legend")
  
  full.dialog <- rk.XML.dialog(
    label = "qIA Analysis",
    rk.XML.tabbook(tabs = list("Basic settings" = list(basic.settings), 
                               "Options Pre-processing" = list(smooth.frame,
                                                               bg.outliers.chk,
                                                               median.chk,
                                                               normalization.frame,
                                                               amptest.chk,
                                                               bg.range.frame,
                                                               threshold.frame, 
                                                               scale.axes.drop), 
                               "Options Quantification" = list(quantification.frame, scale.time.spin),
                               "Plot options" = list(label.frame, generic.plot.options, legend.frame)
    )
    )
  )
  
  JS.calc <- rk.paste.JS(
    echo("options( warn = ", suppress.warnings.chk," )\n"),
    js.selected.x <- rk.JS.vars(selected.x, join = ", "),
    js.var.data <- rk.JS.vars(var.data, join = ", "),
    echo("y <- as.matrix(data.frame(rk.list(", js.var.data,")))\n"),
    ite(single.x.chk, echo("x <- as.matrix(data.frame(rep(rk.list(", js.selected.x,"), ncol(y)))) * ", scale.time.spin,"\n"), 
        echo("x <- as.matrix(data.frame(rk.list(", js.selected.x,"))) * ", scale.time.spin,"\n")),
    echo("if (ncol(x) != ncol(y)) {stop(\"Unequal number of X and Y variables given\")}\n\n"),
    echo(scale.axes.drop),
    # find range of X/Y values needed
    echo("x.range <- range (c (x), na.rm = TRUE)\n"),
    echo("y.range <- range (c (y), na.rm = TRUE)\n\n"),
    
    echo("res.CPP <- lapply(1L:ncol(x), function(i) {\n"),
    echo("if(mean(amptester(y[, i])) > 1) {\n"),
    echo("y.tmp <- log(y[, i])\n"),
    echo("\tres <- CPP(x[, i], y.tmp"),
    ite(smoother.chk, echo(", smoother = FALSE")),
    echo(", method = \"", method.drop,"\""),
    ite(trans.chk, echo(", trans = TRUE, bg.range = c(", bg.range.start.spin,", ", bg.range.end.spin,")")),
    echo(",\n\t\t method.reg = \"", method.reg.drop,"\""),
    ite(bg.outliers.chk, echo(", bg.outliers = TRUE")),
    ite(median.chk, echo(", median = TRUE")),
    echo(",\n\t\t method.norm = \"", method.norm.drop,"\""),
    echo(", qnL = ", qnL.spin,""),
    ite(amptest.chk, echo(",\n amptest = TRUE")),
    ite(manual.chk, echo(", manual = TRUE, nl = ", nL.spin,"")),
    echo("\n)"),
    echo("\t}"),
    echo("})\n"),
    echo("CPP.range <- range(unlist(lapply(1L:ncol(x), function(i) {res.CPP[[i]]$y.norm})), na.rm = TRUE)\n"),
    
    echo("res.cyc.th <- lapply(1L:length(res.CPP), function(i) {\nres <- th.cyc(x[, i], res.CPP[[i]]$y.norm, r = ", r.spin,""),
    ite(auto.chk, echo(", auto = TRUE")),
    ite(linear.chk, echo(", linear = FALSE")),
    echo("\n)})\n\n"),
    echo("res.Cq <- data.frame(Sample = colnames(y),\nCq = sapply(1L:length(res.cyc.th), function (i) {res.cyc.th[[i]][1]}), \nfluo = sapply(1L:length(res.cyc.th), function (i) {res.cyc.th[[i]][2]}))\n"),
    ite(amptest.chk, echo("res.amptest <- t(sapply(1L:length(res.cyc.th), function (i) {data.frame(t(res.CPP[[i]]$y.norm@decisions))}))"))
  )
  
  JS.print <- rk.paste.JS(
    rk.paste.JS.graph(
      echo("type <- rep(\"", pointtype.drop,"\", length.out = length(y))\n"),
      echo("col <- c(1:length(y))\n"),
      echo("pch <- rep(", pch.spin,", length.out = length(y))\n"),
      echo("par(mfrow = c(1,2))\n"),
      echo("plot(NA, NA, xlim = x.range, ylim = y.range, main = \"", plot.main,"\", xlab = \"", plot.xlab,"\", ylab = \"", plot.ylab,"\")\n"),
      # plot variables one X/Y pair at a time
      echo("for (i in 1L:ncol(x)) {\n"),
      echo("\tpoints (x[, i], y[, i], pch = pch[i], type = type[i], col = col[i])\n"),
      echo("\t}\n\n"),
      echo("legend(\"", legend.pos.drop,"\", colnames(y), ncol = ", ncol.legend.spin,", pch = pch, col = col)\n"),
      echo("plot(NA, NA, xlim = x.range, ylim = CPP.range, main = \"", plot.main,"\", xlab = \"", plot.xlab,"\", ylab = \"", plot.ylab,"\")\n"),
      # plot variables one X/Y pair at a time
      echo("lapply(1L:ncol(x), function(i) {points(x[, i], res.CPP[[i]]$y.norm, pch = pch[i], type = type[i], col = col[i])\n"),
      ite(show.th.cyc.chk, echo("\tpoints(as.data.frame(res.cyc.th[[i]]), pch = 4, col = 2)\n")),
      echo("})\n"),
      ite(trans.chk, echo("abline(v = c(", bg.range.start.spin,", ", bg.range.end.spin,"), lty = 2, col = \"grey\")\n")),
      echo("legend(\"", legend.pos.drop,"\", colnames(y), ncol = ", ncol.legend.spin,", pch = pch, col = col)\n"),
      echo("par(mfrow = c(1,1))\n")
    ),
    ite("full", rk.paste.JS(
      echo("\nrk.print(res.Cq)\n"),
      level = 3),
    )
  )
  
  qIAanalysis <<-  rk.plugin.skeleton(
    about = about.info,
    dependencies = dependencies.info,
    xml = list(dialog = full.dialog),
    js = list(require = "chipPCR",
              calculate = JS.calc,
              doPrintout = JS.print,
              results.header = FALSE), # results.header = FALSE is used for backward compatibility (RKWard v. < 0.6.2). Shoudl be removed for RKWard 0.6.3 or later
    rkh = list(plugin.summary, plugin.usage),
    pluginmap = list(
      name = "Quantitative Isothermal Amplification",
      hierarchy = list("analysis", "Quantitative Isothermal Amplification")),
    create=c("pmap","xml","js","desc", "rkh"),
    load = TRUE,
    overwrite = TRUE,
    show = TRUE
  )
})

rk.build.plugin(qIAanalysis)
