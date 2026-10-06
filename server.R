### having example data

# the example file was orginally sequenced from Moriarity Lab work,
  # original file was 3a_PVT1For
examplefile <- "example.ab1"
exampleguide <- "CACTGGAATGACACACGCCC"

# sangerseqR::chromatogram() uses base graphics and restores the graphics parameters
# when it finishes, so the highlight is drawn as a transparent overlay that reproduces
# the geometry of the (single row) chromatogram: same margins, x-limits and trace offsets.
highlightChromatogramBase <- function(obj, trim5, trim3, target.index,
                                      col = "darkblue") {
  avepos <- rowMeans(obj@peakPosMatrix, na.rm = TRUE)
  avepos <- avepos[(1 + trim5):(length(avepos) - trim3)]
  
  last <- min(max(avepos, na.rm = TRUE) + 10, nrow(obj@traceMatrix))
  offset <- max(c(1, avepos[1] - 10))
  tracewidth <- floor(floor(last) - offset) + 1
  rel.pos <- avepos - (offset - 1)
  
  midp <- diff(c(0, rel.pos, tracewidth)) / 2
  k <- target.index - trim5
  
  par(new = TRUE, mar = c(2, 2, 2, 1))
  plot(NA, xlim = c(1, tracewidth), ylim = c(0, 1), type = "n",
       axes = FALSE, xlab = "", ylab = "")
  rect(rel.pos[k] - midp[k], 0, rel.pos[k] + midp[k + 1], 1,
       col = NA, border = col, lwd = 3)
}

# Draws the protospacer chromatogram (top) and the tile plot (bottom) on the open device.
# chromatogram() resets par(fig), so the bottom part of the device is reserved through the
# outer margin and the ggplot is printed into it with a grid viewport.
drawProtospacerFigure <- function(obj, trim5, trim3, tile.plot, target.index = NULL,
                                  fig.height = 8, tile.fraction = 0.4) {
  par(omi = c(fig.height * tile.fraction, 0, 0, 0))
  chromatogram(obj = obj, showcalls = "none", showhets = FALSE,
               trim5 = trim5, trim3 = trim3,
               width = length(obj@primarySeq) - trim5 - trim3)
  if(!is.null(target.index)){
    highlightChromatogramBase(obj, trim5 = trim5, trim3 = trim3, target.index = target.index)
  }
  print(tile.plot, vp = grid::viewport(y = 0, height = tile.fraction, just = "bottom"))
}

shinyServer(
  
  function(input, output) {
    
    
    # reading in file
    input.seqReactive <- reactive({
      # if else things to handle example data loading
      if(input$example) {
        return(readsangerseq(examplefile))
      } else if(!is.null(input$file)) {
        return(readsangerseq(input$file$datapath))
      } else return(validate(
        need(input$file, "Please upload your sanger sequence file")
      ))
      
      
    })
    
    # making basecalls
    
    input.basecallsReactive <- reactive({
      makeBaseCalls(input.seqReactive())
    })
    
    # getting the guide sequence
    guideReactive <- reactive({
      
      
      if(input$example) {
        guide <- (DNAString(exampleguide))
        }
      else{ 
        validate(
          need(input$guide != "", "Please enter a guide RNA sequence")
        )
        
        guide <- (DNAString(input$guide))
        
        if(input$guide.is.reverseComplement){
          guide <- reverseComplement(guide)
        }
      } 
      
      return(guide)
    })
    
    ## pvalue cutoff
    
    p.val.Reactive <- reactive({
      validate(
        need(input$pvalcutoff != "", "Please enter a p-value cutoff")
      )
      
      input$pvalcutoff
    })
    
    
    
    sangsReactive <- reactive({
      input.seq <- input.seqReactive()
      input.basecalls <- input.basecallsReactive()
      input.peakampmatrix <- peakAmpMatrix(input.basecalls)
      
      
      ### creating a sanger sequencing data.frame
      sangs <- CreateSangs(input.peakampmatrix, input.basecalls)
    
      return(sangs)
    })

    
    sangs.filtReactive <- reactive({
      sangs <- sangsReactive()
      
      if(is.na(input$trim3) & is.na(input$trim5)){
        # Remove 5' poor sequencing due to poor primer binding
        sangs.filt <- sangs %>% filter(index > 20)
        # removing crappy end
        peakTotAreaCutoff <- mean(sangs.filt$Tot.area)/10
        sangs.filt %<>% filter(Tot.area > peakTotAreaCutoff)
        return(sangs.filt)
      }else{
        validate(need(!is.na(input$trim3), label = "A value to trim the 3' end "))
        validate(need(!is.na(input$trim5), label = "A value to trim the 5' end "))

        sangs.filt <- sangs[input$trim5:input$trim3, ]
        return(sangs.filt)
      }
      

    })
    
    
    output$trimmedrange <- renderText({
      
      if(!is.na(input$trim5) & !is.na(input$trim3)){
        paste("Trimmed input, starting at ", input$trim5, " and ending at ", input$trim3)
      }

    })
    

    # finding the guide coordinates
    guide.coordReactive <- reactive({
      # the guide coordinates are relative to the index of the sequence that came out of the 
      # makeBaseCalls function
      # the guide is matched to the sequence of the filtered sanger sequencing, so that it doesn't
      # match to any of the low quality regions. 
      
      sangs.filt <- sangs.filtReactive()
      
      # getting the sequence to match to
      filt.sequence <- sangs.filt$base.call %>% paste(collapse = "") %>% DNAString()
      
      # this function finds where the guide matches 
      guide.match <- GetGuideMatch(guideReactive(), filt.sequence)
      
      # Finding the index values
      guide.coord <- list(start = sangs.filt[guide.match$start, "index"],
                          end = sangs.filt[guide.match$end, "index"])
      
      
      return(guide.coord)
      
    })
    
    # index of the base in the sequencing read where editing is expected (NULL if not set/invalid)
    expected.edit.indexReactive <- reactive({
      pos <- input$expected.edit.position
      if(is.null(pos) || length(pos) != 1 || is.na(pos)) return(NULL)
      
      guide.coord <- guide.coordReactive()
      guide.length <- length(guideReactive())
      
      if(pos < 1 || pos > guide.length || pos != round(pos)) return(NULL)
      
      guide.coord$start + as.integer(pos) - 1L
    })
    
    


  
    #######################
    ## Model fitting work
    
    
    nullparams.Reactive <- reactive({
      sangs.filt <- sangs.filtReactive()
      guide.coord <- guide.coordReactive()
      # getting params for the different null models
      null.m.params <- GetNullDistModel(sangs.filt, guide.coord)
    })
    
    
    # Editing reactive
    # this produces a dataframe of the guide region with the probabilities that the nonprimary base shows evidence of editing
    
    editing.Reactive <- reactive({
      # grabbing the reactive stuff
      sangs.filt <- sangs.filtReactive()
      sangs <- sangsReactive()
      guide.coord <- guide.coordReactive()
      guide <- guideReactive()
      null.m.params <- nullparams.Reactive()
      
      
      editing.df <- CreateEditingDF(guide.coord, guide, sangs, null.m.params)
      
      return(editing.df)
    })
    
    
    ### Making QC plots
    
    output$prefilter.totalarea <- renderPlot({
      guide.coord <- guide.coordReactive()
      sangs <- sangsReactive()
      ggplot(sangs, aes(x = index, y = Tot.area)) +
        geom_rect(xmin = guide.coord$start, xmax = guide.coord$end, ymin = 0, ymax = Inf, fill = "lightgrey") +
        geom_line() +
        labs(x = "Base position",
             y = "Total peak area at base position (A + T + C + G)",
             title = "Unfiltered data, total peak area")
    })
    
    output$postfilter.signal.noise <- renderPlotly({
      guide.coord <- guide.coordReactive()
      sangs.filt <- sangs.filtReactive()
      sangs <- sangsReactive()
      
      ## doing some data rearrangement for plotting
      
      catch <- sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
      catch %<>% gather( key = base, value = value,
                         A.area:T.area, A.perc:T.perc) %>%
        separate(col = base, into = c("base", "measure"))
      
      noise <- catch %>% 
        group_by(index, measure) %>% 
        filter(base != base.call) %>%
        summarize(noise = sum(value))
      signal <- catch %>%
        group_by(index, measure) %>%
        filter(base == base.call) %>% 
        summarize(signal = sum(value))
      
      summary.df <- left_join(noise,signal) %>%
        gather(key = type, value = value, noise, signal)
      
      sangs.plot <- sangs.filt %>% dplyr::select(index, base.call) 
      sangs.plot %<>% left_join(summary.df)
      
      ## now acutally making the plot
      
      sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
      
      max.height <- filter(sangs.plot, measure == "area") %>% select(value) %>% max
      
      guide.region <- data.frame(index = c(guide.coord$start, guide.coord$end,
                                           guide.coord$end, guide.coord$start), 
                                 value = c(0, 0, max.height, max.height))
      
      p <- sangs.plot %>% filter(measure == "area") %>%
        ggplot(aes(x = index, y = value)) +
        geom_polygon(data = guide.region, fill = "lightgrey") +
        geom_line(aes(color = type)) +
        scale_color_manual(values = c("#5e3c99", "#e66101")) +  
        labs(title = "Filtered data signal and noise total area",
             x = "Position",
             y = "Peak Area") + 
        theme(legend.title = element_blank())
      ggplotly(p)
      
    })
    
    output$postfilter.noise.perc <- renderPlot({
      guide.coord <- guide.coordReactive()
      sangs.filt <- sangs.filtReactive()
      
      guide.coord <- guide.coordReactive()
      sangs.filt <- sangs.filtReactive()
      
      ## doing some data rearrangement for plotting
      
      catch <- sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
      catch %<>% gather( key = base, value = value,
                         A.area:T.area, A.perc:T.perc) %>%
        separate(col = base, into = c("base", "measure"))
      
      # splitting the catch dataframe into either signal and noise, and calculating the 
      # total noise area or total noise percent
      noise <- catch %>% 
        group_by(index, measure) %>% 
        filter(base != base.call) %>%
        summarize(noise = sum(value))
      signal <- catch %>%
        group_by(index, measure) %>%
        filter(base == base.call) %>% 
        summarize(signal = sum(value))
      
      signal.noise.df <- left_join(noise,signal) %>%
        gather(key = type, value = value, noise, signal)
      
      # making the plotting df
      sangs.plot <- sangs.filt %>% dplyr::select(index, base.call) 
      sangs.plot %<>% left_join(signal.noise.df)
      sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
      
      sangs.plot %>% filter(measure == "perc", type == "noise") %>%
        ggplot(aes(x = index, y = value)) +
        geom_area(aes(fill = type)) +
        scale_fill_manual(values = c("#e66101")) + 
        annotate("rect", xmin=guide.coord$start, xmax=guide.coord$end,
                 ymin=0, ymax=Inf, alpha=1/5, fill="black") +
        guides(fill = FALSE) + 
        labs(title = "Percent peak area noise",
             x = "Position",
             y = "Percent peak area noise")
    })
    
    
    output$chromatogram <- renderPlot({
      guide.coord <- guide.coordReactive()
      input.basecalls <- input.basecallsReactive()
      
      chromatogram(obj = input.basecalls, 
                   trim5 = (guide.coord$start-1), trim3 = length(input.basecalls@primarySeq) - guide.coord$end,
                   width = (guide.coord$end - guide.coord$start + 1))
      
      target.index <- expected.edit.indexReactive()
      if(!is.null(target.index)){
        highlightChromatogramBase(input.basecalls, trim5 = (guide.coord$start-1),
                                  trim3 = length(input.basecalls@primarySeq) - guide.coord$end,
                                  target.index = target.index)
      }
    })
    
    
    base.infoReactive <- reactive({
      # this reactive is collecting information on each base:
      # - average percent signal
      # critical value for each base for the cutoff of significance
      
      sangs.filt <- sangs.filtReactive()
      null.m.params <- nullparams.Reactive()
      p.val.cutoff <- p.val.Reactive()
      
      # finding the average percent signal for each base
      avg.base <- sangs.filt %>% gather(key = focal.base, value = value, 
                                A.area:T.area, A.perc:T.perc) %>%
        separate(col = focal.base, into = c("focal.base", "measure")) %>% 
        spread(key = measure, value = value) %>% 
        filter(base.call == focal.base) %>% 
        group_by(focal.base) %>% 
        summarize(avg.percsignal = mean(perc),
                  avg.areasignal = mean(area))
      
      # getting the critical value for each base
      crit.vals <- c(
        a = qZAGA(p = p.val.cutoff, mu = null.m.params$a$mu,
            sigma = null.m.params$a$sigma,
            nu = null.m.params$a$nu,
            lower.tail = FALSE),
        c = qZAGA(p = p.val.cutoff, mu = null.m.params$c$mu,
                  sigma = null.m.params$c$sigma,
                  nu = null.m.params$c$nu,
                  lower.tail = FALSE),
        g = qZAGA(p = p.val.cutoff, mu = null.m.params$g$mu,
                  sigma = null.m.params$g$sigma,
                  nu = null.m.params$g$nu,
                  lower.tail = FALSE),
        t = qZAGA(p = p.val.cutoff, mu = null.m.params$t$mu,
                  sigma = null.m.params$t$sigma,
                  nu = null.m.params$t$nu,
                  lower.tail = FALSE)
        )
      # getting fillibens
      fil <- lapply(null.m.params, FUN = function(x){x$fillibens})
      filvec <- c(a = fil$a, c = fil$c, g = fil$g, t = fil$t)
      
      # getting mu -- a measure of dispersion
      mul <- lapply(null.m.params, FUN = function(x){x$mu})
      mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
      
      return(data.frame(avg.base, crit.perc.area = crit.vals, mu = mulvec, fillibens = filvec))
    })
    
    
    output$baseinfo.table <- renderTable({
      base.info <- base.infoReactive()
      
      temp <- base.info
      names(temp) <- c("Base", "Average percent signal",
                       "Average peak area",  "Critical percent value",
                       "model mu",  "Fillibens correlation")
      row.names(temp) <- NULL
      
      return(temp)
    })
    
    ##############################################
    output$chromatogram_two <- renderPlot({
      guide.coord <- guide.coordReactive()
      input.basecalls <- input.basecallsReactive()
      chromatogram(obj = input.basecalls,
                   showcalls = "none",
                   showhets = FALSE,
                   trim5 = (guide.coord$start-1), trim3 = length(input.basecalls@primarySeq) - guide.coord$end,
                   width = (guide.coord$end - guide.coord$start + 1))
      
      target.index <- expected.edit.indexReactive()
      if(!is.null(target.index)){
        highlightChromatogramBase(input.basecalls, trim5 = (guide.coord$start-1),
                                  trim3 = length(input.basecalls@primarySeq) - guide.coord$end,
                                  target.index = target.index)
      }
    })
    
    output$editing.table.plot <- renderPlot({
      editing.table.ggReactive()
    })
    
    editing.table.ggReactive <- reactive({
      editing.df <- editing.Reactive()
      sangs.filt <- sangs.filtReactive()
      null.m.params <- nullparams.Reactive()
      p.val.cutoff <- p.val.Reactive()
      
      #### Repeat code for getting avg.base from base.infoReactive
      # finding the average percent signal for each base
      avg.base <- sangs.filt %>% gather(key = focal.base, value = value, 
                                        A.area:T.area, A.perc:T.perc) %>%
        separate(col = focal.base, into = c("focal.base", "measure")) %>% 
        spread(key = measure, value = value) %>% 
        filter(base.call == focal.base) %>% 
        group_by(focal.base) %>% 
        summarize(avg.percsignal = mean(perc),
                  avg.areasignal = mean(area))
      
      # finding the model mu
      mul <- lapply(null.m.params, FUN = function(x){x$mu})
      mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
      ####
      
      ### Reshape data
      
      edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                                         A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
        separate(col = focal.base, into = c("focal.base", "measure"))
      
      edit.spread <- edit.long %>% 
        spread(key = measure, value = value) 
      
      color.cutoff = min(avg.base$avg.percsignal - mulvec)
      edit.color <- edit.spread %>% 
        mutate(adj.perc = {ifelse(perc >= color.cutoff,
                                  100,
                                  perc)
        } %>% as.numeric) %>%
        filter(pval < p.val.cutoff)
      

      
      # unfilled box around the column of the expected editing position
      target.index <- expected.edit.indexReactive()
      box.layer <- NULL
      if(!is.null(target.index)){
        col.pos <- match(target.index, sort(unique(edit.spread$index)))
        if(!is.na(col.pos)){
          box.layer <- annotate("rect", xmin = col.pos - 0.5, xmax = col.pos + 0.5,
                                ymin = 0.5, ymax = 4.5,
                                fill = NA, colour = "darkblue", linewidth = 1.2)
        }
      }
      
      #### make editing_table
      if(any(edit.color$adj.perc != 100)){
edit.spread %>%
        ggplot(aes(x = as.factor(index), y = focal.base)) + 
        geom_tile(data = edit.color, aes(fill = adj.perc)) + 
        geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
        box.layer +
        guides(fill = FALSE) + 
        scale_fill_continuous(low = "#f7a8a8", high = "#9acdee") +
        scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
        labs(x = NULL, y = NULL) + 
        theme(axis.ticks = element_blank(),
              axis.text=element_text(size=16),
              plot.title = element_text(hjust = 0, size = 16),
              plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
              panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
              plot.background = element_rect(fill = "transparent",colour = NA)
        ) +
        coord_fixed(1)} 
      else
      {edit.spread %>%
        ggplot(aes(x = as.factor(index), y = focal.base)) + 
        geom_tile(data = edit.color, fill = "#9acdee") + 
        geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
        box.layer +
        guides(fill = FALSE) + 
        scale_x_discrete(position = "top", labels = editing.df$guide.seq) +
        labs(x = NULL, y = NULL) + 
        theme(axis.ticks = element_blank(),
              axis.text=element_text(size=16),
              plot.title = element_text(hjust = 0, size = 16),
              plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
              panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
              plot.background = element_rect(fill = "transparent",colour = NA)
        ) +
        coord_fixed(1)}
    })
    
    # export of the protospacer chromatogram together with the percent-area tile plot
    protospacerFigureDevice <- function(file, type){
      guide.coord <- guide.coordReactive()
      input.basecalls <- input.basecallsReactive()
      tile.plot <- editing.table.ggReactive()
      target.index <- expected.edit.indexReactive()
      
      fig.width <- 14
      fig.height <- 8
      if(type == "png"){
        png(file, width = fig.width, height = fig.height, units = "in", res = 200, bg = "white")
      } else {
        svg(file, width = fig.width, height = fig.height, bg = "white")
      }
      on.exit(dev.off())
      
      drawProtospacerFigure(obj = input.basecalls,
                            trim5 = (guide.coord$start-1),
                            trim3 = length(input.basecalls@primarySeq) - guide.coord$end,
                            tile.plot = tile.plot,
                            target.index = target.index,
                            fig.height = fig.height)
    }
    
    output$downloadProtospacerPNG <- downloadHandler(
      filename = function(){ "EditR_protospacer.png" },
      content = function(file){ protospacerFigureDevice(file, "png") }
    )
    
    output$downloadProtospacerSVG <- downloadHandler(
      filename = function(){ "EditR_protospacer.svg" },
      content = function(file){ protospacerFigureDevice(file, "svg") }
    )
    
###########################################
    output$editing.quad.plot <- renderPlot({
      editing.df <- editing.Reactive()
      null.m.params <- nullparams.Reactive()
      p.val.cutoff <- p.val.Reactive()
      
      edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                                         A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
        separate(col = focal.base, into = c("focal.base", "measure"))
      
      p.a <- makeEditingBarPlot(edit.long = edit.long, null.m.params = null.m.params$a,
                                base = "A", pval = p.val.cutoff, editing.df)
      p.c <- makeEditingBarPlot(edit.long, null.m.params$c,
                                base = "C", pval = p.val.cutoff, editing.df)
      p.g <- makeEditingBarPlot(edit.long, null.m.params$g,
                                base = "G", pval = p.val.cutoff, editing.df)
      p.t <- makeEditingBarPlot(edit.long, null.m.params$t,
                                base = "T", pval = p.val.cutoff, editing.df)
      
      grid.arrange(p.a, p.c, p.g, p.t)
    })



##########################33
    output$editingtable <- renderTable(digits = -3, expr = { # maybe need to use renderMarkdown?
      editing.df <- editing.Reactive()
      p.val.cutoff <- p.val.Reactive()
      
      edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                                         A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
        separate(col = focal.base, into = c("focal.base", "measure"))
      
      edit.spread <- edit.long %>% spread(key = measure, value = value)
    
      
      editingtable <- edit.spread %>% arrange(guide.position) %>% select(index,guide.position, guide.seq, base.call, focal.base, perc, pval) %>% 
        mutate(index = as.character(index), 
               guide.position = as.character(guide.position),
               perc = format(perc, digits = 2),
               signficant = ifelse(pval < p.val.cutoff, yes = "*", no = ""))
      
      names(editingtable) <-  c("Sanger position", "Guide position","Guide sequence", "Sanger base call",
                         "Focal base", "Focal base peak area", "p value", "")
      
      return(editingtable)
    })
    
    output$modelfits <- renderTable({
      null.m.params <- nullparams.Reactive()
      
      temp <- data.frame(`Base` = c("A", "C", "G", "T"),
                 `Fillibens Correlation` = c(null.m.params$a$fillibens, 
                                          null.m.params$g$fillibens, 
                                          null.m.params$c$fillibens, 
                                          null.m.params$t$fillibens))
      
    names(temp) <- c("Base", "Fillibens Correlation")
    return(temp)
    })
    
    
#################################3
    
    output$downloadReport <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "report.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
         # pulling in the reactive things to generate the report with:
        
        params <- list(null.m.params = nullparams.Reactive(),
                       guide = guideReactive(),
                       guide.coord = guide.coordReactive(),
                       editing.df = editing.Reactive(),
                       p.val.cutoff = p.val.Reactive(),
                       input.seq = input.seqReactive(),
                       sangs = sangsReactive(),
                       sangs.filt = sangs.filtReactive(),
                       base.info = base.infoReactive(),
                       editrversion = editrversion)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )

})
