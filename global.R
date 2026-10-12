# version
editrversion = "1.0.10"

# loading all the packages into the global environment

library(shiny)
library(Biostrings)
library(sangerseqR)
library(gamlss)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(rmarkdown)
library(plotly)
library(yaml)


## Defining functions



GetGuideMatch <- function(guide, input.seq)
  # this function takes in the guide sequence, and the primary seq from a makeBaseCalls object
  # guide should be a DNAString object -- DNAString(guideinput)
  # input.seq = input.basecalls@primarySeq
  # returns guide.coord, a list with $match (forward / reverse), $start, and $end
  # and finds where the guide matches
  # it looks for both a forward and reverse match
  # in deciding if it's forward / reverse, it takes the match that is the longest
  # This function tries to extend the guide coordinates in case there is a partial match -- this might break
  # match to guide
  {

  
  # alignments
  align.f  <- pairwiseAlignment(pattern=guide, subject=input.seq, type="overlap")

  guide.coord <- list(match = "forward", start = align.f@subject@range@start,
                      end = align.f@subject@range@start + align.f@subject@range@width-1)
  
  # check to see if a match was actually found
  if(align.f@pattern@range@width <= 1){
    stop("Failed to find a forward match of the gRNA -- should it be a reverse complement?")
    return(NA)
  }
  

  # check to see if the matched region == the guide length
  is.match.guide.length <- (guide.coord$end - guide.coord$start + 1) == length(guide)
  
  # otherwise figure out what part of the guide matched to where, and update the guide region to include the unmatched areas. 
  if(is.match.guide.length != TRUE){
    
    # find the diff from beginning
    startdiff <- 1 - align.f@pattern@range@start
    # find the diff from the end
    enddiff <- 20 - (align.f@pattern@range@start + align.f@pattern@range@width-1)
    # add these numbers to the guide.cord
    guide.coord$start <- guide.coord$start + startdiff
    guide.coord$end <- guide.coord$end + enddiff
    
    # check to make sure that I didn't do something bad
    # if(guide.coord$start < 1 | guide.coord$end > length(guide)){
    #   stop("Something went wrong in matching the guide to the sequence, partial match of the guide found.
    #          Error in extending the match to the full guide length.")
    # }
    
    
  }
  
  return(guide.coord)
}


CreateSangs <- function(peakAmp, basecalls){
  # this function takes the peak amp matrix, and the input basecalls, and creates an 
  # unfiltered sanger sequencing dataframe (that I like to call sangs)
  
  
  sangs <- as.data.frame(peakAmp)
  names(sangs) <- c("A.area","C.area","G.area","T.area")
  
  # adding total area, and then calculating the percent area of each base
  sangs %<>% 
    mutate(Tot.area = A.area + C.area + G.area + T.area,
           A.perc = 100*A.area / Tot.area,
           C.perc = 100*C.area / Tot.area,
           G.perc = 100*G.area / Tot.area,
           T.perc = 100*T.area / Tot.area) 
  
  # adding on base calls
  sangs$base.call <- strsplit(x = toString(basecalls@primarySeq), split = "") %>%
    unlist
  
  # adding an index
  sangs$index <- seq_along(sangs$base.call)
  return(sangs)
}


GetNullDistModel <- function(sangs.filt, guide.coord)
  # this function takes the sangs.filt dataframe, gathers values for the null distribution
  # for each base, and then fits a zero-adjusted gamma distribution to these values
  # it returns a list of data.frames, each with the parameters of a null distribution 
  # for each base
  {

  
  sangs.filt %<>% filter(!(index %in% (guide.coord$start:guide.coord$end)) )
  
  nvals <- list()
  nvals$t <- sangs.filt %>% filter(base.call != "T") %>% select(T.perc) %>% unlist()
  nvals$c <- sangs.filt %>% filter(base.call != "C") %>% select(C.perc) %>% unlist()
  nvals$g <- sangs.filt %>% filter(base.call != "G") %>% select(G.perc) %>% unlist()
  nvals$a <- sangs.filt %>% filter(base.call != "A") %>% select(A.perc) %>% unlist()
  
  # Updated 3.26.19 to account for ultra clean sequencing
  replacement_zaga = c(rep(0, 989), 0.00998720389310502, 0.00998813447664401,0.009992887520785,
                       0.00999585366068316, 0.00999623914632598, 0.00999799013526835, 0.010001499423723,
                       0.0100030237039207, 0.0100045782875701, 0.0100048452355807, 0.0100049548867042)
  
  n.models <- lapply(nvals, FUN = function(x){
    set.seed(1)
    if((unique(x)[1] == 0 & length(unique(x)) == 1) |
       (unique(x)[1] == 0 & length(unique(x)) == 2 & table(x)[2] == 1))
    {x = replacement_zaga; message("Replacement vector used for low noise.")} # add noise if all 0s, or all 0s and one other value.
    tryCatch(gamlss((x)~1, family = ZAGA), error=function(e) # Progressively step up the mu.start if it fails
      tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 1), error=function(e) 
        tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 2), error=function(e) 
          tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 3), error=function(e) # additional step added.
            gamlss((x)~1, family = ZAGA, mu.start = mean(x))
          )
        )
      )
    )
    # throws errors when a completely 0 vector
  })
  
  null.m.params <- lapply(n.models, FUN = function(x){
    mu <- exp(x$mu.coefficients[[1]])
    sigma <- exp(x$sigma.coefficients[[1]])
    nu.logit <- x$nu.coefficients[[1]]
    nu <- exp(nu.logit)/(1+exp(nu.logit))
    fillibens <-cor(as.data.frame(qqnorm(x$residuals, plot = FALSE)))[1,2]
    
    return(data.frame(mu= mu, sigma = sigma, nu = nu, fillibens = fillibens))
  })
  
  return(null.m.params)
}


CreateEditingDF <- function(guide.coord, guide, sangs, null.m.params){
  # this function creates the editing.df that contains the prob that the peak areas 
  # are not part of the noise distibution
  
  guide.df <- sangs[guide.coord$start:guide.coord$end,]
  
  # adding the guide sequence onto the data.frame
  guide.df$guide.seq <- guide %>% toString() %>% strsplit(. , "") %>% unlist
  
  
  # usage (params = null.m.params$t, perc = guide.df$T.perc)
  calcBaseProb <- function(params, perc){
    outprobs <- pZAGA(q = perc,
                      mu = params$mu,
                      sigma = params$sigma,
                      nu = params$nu,
                      lower.tail = FALSE)
  }
  
  guide.df$T.pval <- calcBaseProb(params = null.m.params$t, perc = guide.df$T.perc)
  guide.df$C.pval <- calcBaseProb(params = null.m.params$c, perc = guide.df$C.perc)
  guide.df$G.pval <- calcBaseProb(params = null.m.params$g, perc = guide.df$G.perc)
  guide.df$A.pval <- calcBaseProb(params = null.m.params$a, perc = guide.df$A.perc)
  
  guide.df$guide.position <- seq_along(guide.df$guide.seq)
  
  editing.df <- guide.df  
  return(editing.df)
}


makeEditingBarPlot <- function(edit.long, null.m.params, base, pval, editing.df){
  # makeEditingBarPlot(edit.filt = filter(editing.df, guide.seq != "T"),
  # null.m.params = null.m.params$t, base = "T, pval = p.val.cutoff)
  
  
  cutoff.line <- qZAGA(p = pval, mu = null.m.params$mu,
                       sigma = null.m.params$sigma, 
                       nu = null.m.params$nu, lower.tail = FALSE)
  
  edit.filt <- edit.long %>% filter(guide.seq != base, focal.base == base)
  
  edit.filt %>% filter(measure == "perc") %>% 
    ggplot(aes(x = guide.position, y = value)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    coord_cartesian(ylim = c(0,100)) + 
    scale_x_continuous(breaks = seq_along(editing.df$guide.seq)) + 
    annotate(geom = "segment", x = 0, xend = length(editing.df$guide.seq),
             y = cutoff.line, yend = cutoff.line) +
    labs(title = paste0("Percent ", base, " editing"),
         x = "Base postion in guide RNA",
         y = "Percent Total peak area") + 
    theme(axis.text.x = element_text(size = 8))
  
}



GUIDE_NOT_FOUND <- "The guide sequence was not found in the sequencing read. Check the guide, the 'Guide sequence is reverse complement' box and the trim values."

AnalyseSangerSample <- function(seq, guide, trim5, trim3){
  # runs the same analysis as for the edited sample on any sangerseq object
  # (used for every uploaded file); returns everything needed for plotting.
  # Problems are reported with stop() so that a bad file does not break the other samples.
  basecalls <- makeBaseCalls(seq)
  sangs <- CreateSangs(peakAmpMatrix(basecalls), basecalls)
  
  if(is.na(trim5) & is.na(trim3)){
    sangs.filt <- sangs %>% filter(index > 20)
    peakTotAreaCutoff <- mean(sangs.filt$Tot.area)/10
    sangs.filt %<>% filter(Tot.area > peakTotAreaCutoff)
  } else {
    if(is.na(trim5) || is.na(trim3)) stop("Values to trim both the 5' and 3' end must be provided")
    if(!(trim5 >= 1 && trim3 <= nrow(sangs) && trim5 < trim3)){
      stop(paste0("Trim values must satisfy 1 <= 5' start < 3' end <= ", nrow(sangs),
                  " (the number of called bases in this file)"))
    }
    sangs.filt <- sangs[trim5:trim3, ]
  }
  
  filt.sequence <- sangs.filt$base.call %>% paste(collapse = "") %>% DNAString()
  guide.match <- tryCatch(GetGuideMatch(guide, filt.sequence), error = function(e) NULL)
  if(is.null(guide.match) || length(guide.match$start) != 1 || is.na(guide.match$start) || is.na(guide.match$end)) stop(GUIDE_NOT_FOUND)
    guide.coord <- list(start = sangs.filt[guide.match$start, "index"],
                        end = sangs.filt[guide.match$end, "index"])
    if(anyNA(guide.coord)) stop(GUIDE_NOT_FOUND)
  null.m.params <- GetNullDistModel(sangs.filt, guide.coord)
  editing.df <- CreateEditingDF(guide.coord, guide, sangs, null.m.params)
  
  list(basecalls = basecalls, sangs = sangs, sangs.filt = sangs.filt,
       guide.coord = guide.coord, null.m.params = null.m.params,
       editing.df = editing.df)
}


NormaliseToWildType <- function(edited.df, wt.df, het.position = NULL){
  # Returns long data.frame (guide.position, guide.seq, focal.base, perc) with the
  # normalised percent area of each base in the guide.
  # - at the heterozygous position: Abbott formula, (edited - wt) / (100 - wt) * 100
  # - all other positions: PLACEHOLDER, simple background subtraction max(edited - wt, 0)
  bases <- c("A", "C", "G", "T")
  e <- as.matrix(edited.df[, paste0(bases, ".perc")])
  w <- as.matrix(wt.df[, paste0(bases, ".perc")])
  
  norm <- pmax(e - w, 0)
  if(!is.null(het.position)){
    i <- het.position
    abbott <- (e[i, ] - w[i, ]) / pmax(100 - w[i, ], 1e-6) * 100
    norm[i, ] <- pmin(pmax(abbott, 0), 100)
  }
  colnames(norm) <- bases
  
  data.frame(guide.position = rep(edited.df$guide.position, times = 4),
             guide.seq = rep(edited.df$guide.seq, times = 4),
             focal.base = rep(bases, each = nrow(norm)),
             perc = as.vector(norm),
             stringsAsFactors = FALSE)
}


######################
## Figures: chromatogram + base-area table, shared by wild-type, edited and normalised

# same colours as sangerseqR::chromatogram
TraceColours <- c(A = "green", C = "blue", G = "black", T = "red")

GuideTraceData <- function(basecalls, guide.coord){
  # trace of the guide region, x rescaled so that the peak of guide base k sits at x = k
  traces <- basecalls@traceMatrix
  avepos <- rowMeans(basecalls@peakPosMatrix, na.rm = TRUE)[guide.coord$start:guide.coord$end]
  n <- length(avepos)
  edge5 <- avepos[1] - (avepos[2] - avepos[1]) / 2
  edge3 <- avepos[n] + (avepos[n] - avepos[n - 1]) / 2
  rng <- max(1, ceiling(edge5)):min(nrow(traces), floor(edge3))
  x <- approx(c(edge5, avepos, edge3), c(0.5, seq_len(n), n + 0.5), xout = rng, rule = 2)$y
  
  do.call(rbind, lapply(1:4, function(b){
    data.frame(x = x, y = traces[rng, b], base = c("A", "C", "G", "T")[b])
  }))
}

NormalisedTraceData <- function(norm.long){
  # idealised chromatogram: one gaussian peak per base and position, height = normalised percent area
  grid <- seq(-0.5, 0.5, length.out = 41)
  do.call(rbind, lapply(seq_len(nrow(norm.long)), function(r){
    data.frame(x = norm.long$guide.position[r] + grid,
               y = norm.long$perc[r] * exp(-grid^2 / (2 * 0.15^2)),
               base = norm.long$focal.base[r],
               id = r)
  }))
}

GuideChromatogramPlot <- function(trace.df, n.positions, highlight = NULL,
                                  y.label = NULL, y.max = NULL){
  p <- ggplot(trace.df, aes(x = x, y = y, colour = base,
                            group = if("id" %in% names(trace.df)) id else base))
  if(!is.null(highlight)){
    p <- p + annotate("rect", xmin = highlight - 0.5, xmax = highlight + 0.5,
                      ymin = -Inf, ymax = Inf, fill = NA, colour = "darkblue", linewidth = 1.2)
  }
  p <- p +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = TraceColours, guide = "none") +
    scale_x_continuous(limits = c(0.5, n.positions + 0.5), expand = c(0, 0)) +
    labs(x = NULL, y = y.label) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank())
  
  if(is.null(y.max)){
    p + coord_cartesian(ylim = c(0, max(trace.df$y) * 1.05)) +
      theme(axis.text.y = element_blank())
  } else {
    p + coord_cartesian(ylim = c(0, y.max))
  }
}

EditingTableData <- function(editing.df, sangs.filt, null.m.params, p.val.cutoff){
  avg.base <- sangs.filt %>% gather(key = focal.base, value = value,
                                    A.area:T.area, A.perc:T.perc) %>%
    separate(col = focal.base, into = c("focal.base", "measure")) %>%
    spread(key = measure, value = value) %>%
    filter(base.call == focal.base) %>%
    group_by(focal.base) %>%
    summarize(avg.percsignal = mean(perc),
              avg.areasignal = mean(area))
  
  mul <- lapply(null.m.params, FUN = function(x){x$mu})
  mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
  
  edit.spread <- editing.df %>% gather(key = focal.base, value = value,
                                       A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
    separate(col = focal.base, into = c("focal.base", "measure")) %>%
    spread(key = measure, value = value)
  
  color.cutoff <- min(avg.base$avg.percsignal - mulvec)
  edit.color <- edit.spread %>%
    mutate(adj.perc = as.numeric(ifelse(perc >= color.cutoff, 100, perc))) %>%
    filter(pval < p.val.cutoff)
  
  list(spread = edit.spread, color = edit.color, color.cutoff = color.cutoff,
       guide.seq = editing.df$guide.seq)
}

NormalisedTableData <- function(norm.long, edited.table.data, p.val.cutoff){
  # coloured where the edited sample is significantly different from its noise and the
  # normalised value is above 0
  sig <- edited.table.data$spread %>% dplyr::select(guide.position, focal.base, pval)
  edit.spread <- norm.long %>% left_join(sig, by = c("guide.position", "focal.base"))
  edit.color <- edit.spread %>%
    filter(pval < p.val.cutoff, perc > 0) %>%
    mutate(adj.perc = as.numeric(ifelse(perc >= edited.table.data$color.cutoff, 100, perc)))
  
  list(spread = edit.spread, color = edit.color, guide.seq = edited.table.data$guide.seq)
}

GuideTablePlot <- function(table.data, highlight = NULL){
  positions <- sort(unique(table.data$spread$guide.position))
  p <- ggplot(table.data$spread, aes(x = guide.position, y = focal.base))
  if(any(table.data$color$adj.perc != 100)){
    p <- p + geom_tile(data = table.data$color, aes(fill = adj.perc)) +
      scale_fill_continuous(low = "#f7a8a8", high = "#9acdee", guide = "none")
  } else {
    p <- p + geom_tile(data = table.data$color, fill = "#9acdee")
  }
  if(!is.null(highlight)){
    p <- p + annotate("rect", xmin = highlight - 0.5, xmax = highlight + 0.5, ymin = 0.5, ymax = 4.5,
                      fill = NA, colour = "darkblue", linewidth = 1.2)
  }
  p +
    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +
    scale_x_continuous(position = "top", breaks = positions, labels = table.data$guide.seq,
                       limits = c(0.5, max(positions) + 0.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(size = 14),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA))
}

GuideFigure <- function(chromatogram.plot, table.plot, title){
  # chromatogram on top, table below, panels aligned so that positions line up
  body <- cowplot::plot_grid(chromatogram.plot, table.plot, ncol = 1,
                             align = "v", axis = "lr", rel_heights = c(1.2, 1))
  heading <- cowplot::ggdraw() +
    cowplot::draw_label(title, x = 0.01, hjust = 0, fontface = "bold", size = 16)
  cowplot::plot_grid(heading, body, ncol = 1, rel_heights = c(0.07, 1)) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
}

SaveGuideFigure <- function(figure, file, format = c("png", "svg")){
  format <- match.arg(format)
  if(format == "png"){
    ggsave(file, plot = figure, device = "png", width = 12, height = 6.5, dpi = 300, bg = "white")
  } else {
    ggsave(file, plot = figure, width = 12, height = 6.5,
           device = function(filename, width, height, ...) grDevices::svg(filename, width, height))
  }
}

EDITING_BYSTANDER_THRESHOLD <- 10

EditingSummaryLines <- function(perc.list, guide.seq, sample.names, position, threshold = EDITING_BYSTANDER_THRESHOLD){
  # on-target editing at `position` and bystander editing elsewhere in the guide.
  # perc.list: one matrix per sample (rows = guide positions, columns A, C, G, T) with the percent
  # of each base - raw percent areas, or Abbott / non-treated normalised values.
  # A change is any base that differs from the guide base; with several samples the
  # percentages are averaged across samples.
  # returns list(on.target = chr, off.target = chr)
  bases <- c("A", "C", "G", "T")
  n <- length(perc.list)
  mean.perc <- Reduce(`+`, perc.list) / n
  
  # guide base -> NA
  nonref <- function(m){
    for(i in seq_along(guide.seq)) m[i, bases == toupper(guide.seq[i])] <- NA
    m
  }
  fmt <- function(x) sprintf("%.1f%%", x)
  
  mean.nonref <- nonref(mean.perc)
  ref.base <- toupper(guide.seq[position])
  alt <- mean.nonref[position, ]
  main <- which.max(alt)
  on.target <- paste0("On-target editing at guide position ", position, " (", ref.base, " to ", bases[main], "): ",
                      fmt(alt[main]))
  if(n > 1){
    each <- vapply(perc.list, function(m) nonref(m)[position, main], numeric(1))
    on.target <- paste0(on.target, " (mean of ", n, " samples; ",
                        paste0(sample.names, ": ", fmt(each), collapse = ", "), ")")
  }
  
  mean.nonref[position, ] <- NA
  hits <- which(!is.na(mean.nonref) & mean.nonref >= threshold, arr.ind = TRUE)
  hits <- hits[order(hits[, 1], hits[, 2]), , drop = FALSE]
  off.target <- if(nrow(hits) == 0){
    paste0("No bystander editing of ", threshold, "% or more detected at other guide positions.")
  } else {
    paste0("Bystander editing (", threshold, "% or more) at guide position ", hits[, 1], " (",
           toupper(guide.seq[hits[, 1]]), " to ", bases[hits[, 2]], "): ", fmt(mean.nonref[hits]))
  }
  list(on.target = on.target, off.target = off.target)
}

RawPercMatrix <- function(editing.df){
  as.matrix(editing.df[, paste0(c("A", "C", "G", "T"), ".perc")])
}

NormalisedPercMatrix <- function(edited.df, nt.df, het.position){
  # Abbott-normalised percent of each base (see NormaliseToWildType) as a positions x bases matrix
  norm.long <- NormaliseToWildType(edited.df, nt.df, het.position)
  matrix(norm.long$perc, ncol = 4, dimnames = list(NULL, c("A", "C", "G", "T")))
}
