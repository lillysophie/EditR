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

shinyServer(
  
  function(input, output) {
    
    
    # reading in file
    input.seqReactive <- reactive({
      # if else things to handle example data loading
      files <- edited.filesReactive()
      readsangerseq(files$datapath[selectedSampleReactive()])
      
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
        validate(need(input$trim5 >= 1 && input$trim3 <= nrow(sangs) && input$trim5 < input$trim3,
                      paste0("Trim values must satisfy 1 <= 5' start < 3' end <= ", nrow(sangs),
                             " (the number of called bases)")))
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
      guide.match <- tryCatch(GetGuideMatch(guideReactive(), filt.sequence), error = function(e) NULL)
      validate(need(!is.null(guide.match) && length(guide.match$start) == 1 && !is.na(guide.match$start) && !is.na(guide.match$end),
                    "The guide sequence was not found in the sequencing read. Check the guide, the 'Guide sequence is reverse complement' box and the trim values."))
      
      # Finding the index values
      guide.coord <- list(start = sangs.filt[guide.match$start, "index"],
                          end = sangs.filt[guide.match$end, "index"])
      validate(need(!anyNA(guide.coord),
                    "The guide sequence was not found in the sequencing read. Check the guide, the 'Guide sequence is reverse complement' box and the trim values."))
      
      
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
    
    
    ### Making QC plots (one set per uploaded edited file, built from edited.analysesReactive)
    
    qcPrefilterPlot <- function(a){
      ggplot(a$sangs, aes(x = index, y = Tot.area)) +
        geom_rect(xmin = a$guide.coord$start, xmax = a$guide.coord$end, ymin = 0, ymax = Inf, fill = "lightgrey") +
        geom_line() +
        labs(x = "Base position",
             y = "Total peak area at base position (A + T + C + G)",
             title = "Unfiltered data, total peak area")
    }
    
    qcSignalNoiseLong <- function(a){
      catch <- a$sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
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
      
      summary.df <- left_join(noise, signal) %>%
        gather(key = type, value = value, noise, signal)
      
      sangs.plot <- a$sangs.filt %>% dplyr::select(index, base.call) 
      sangs.plot %<>% left_join(summary.df)
      sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
      sangs.plot
    }
    
    qcSignalNoiseGg <- function(a){
      sangs.plot <- qcSignalNoiseLong(a)
      guide.coord <- a$guide.coord
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
      p
    }
    
    qcSignalNoisePlot <- function(a) ggplotly(qcSignalNoiseGg(a))
    
    qcNoisePercPlot <- function(a){
      guide.coord <- a$guide.coord
      qcSignalNoiseLong(a) %>% filter(measure == "perc", type == "noise") %>%
        ggplot(aes(x = index, y = value)) +
        geom_area(aes(fill = type)) +
        scale_fill_manual(values = c("#e66101")) + 
        annotate("rect", xmin=guide.coord$start, xmax=guide.coord$end,
                 ymin=0, ymax=Inf, alpha=1/5, fill="black") +
        guides(fill = "none") + 
        labs(title = "Percent peak area noise",
             x = "Position",
             y = "Percent peak area noise")
    }
    
    qcChromatogram <- function(a){
      guide.coord <- a$guide.coord
      basecalls <- a$basecalls
      trim3 <- length(basecalls@primarySeq) - guide.coord$end
      chromatogram(obj = basecalls, 
                   trim5 = (guide.coord$start-1), trim3 = trim3,
                   width = (guide.coord$end - guide.coord$start + 1))
      
      pos <- expectedPositionReactive()
      if(!is.null(pos)){
        highlightChromatogramBase(basecalls, trim5 = (guide.coord$start-1), trim3 = trim3,
                                  target.index = guide.coord$start + pos - 1L)
      }
    }
    
    baseInfoFor <- function(a){
      # this reactive is collecting information on each base:
      # - average percent signal
      # critical value for each base for the cutoff of significance
      
      sangs.filt <- a$sangs.filt
      null.m.params <- a$null.m.params
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
    }
    
    base.infoReactive <- reactive({
      baseInfoFor(list(sangs.filt = sangs.filtReactive(), null.m.params = nullparams.Reactive()))
    })
    
    
    baseInfoTable <- function(base.info){
      temp <- base.info
      temp <- base.info
      names(temp) <- c("Base", "Average percent signal",
                       "Average peak area",  "Critical percent value",
                       "model mu",  "Fillibens correlation")
      row.names(temp) <- NULL
      
      return(temp)
    }
    
    ##############################################
    ## Uploaded files: several edited samples, one non-treated control
    
    edited.filesReactive <- reactive({
      if(input$example){
        return(data.frame(name = "example.ab1", datapath = examplefile, stringsAsFactors = FALSE))
      }
      validate(need(!is.null(input$file), "Please upload your sanger sequence file(s)"))
      data.frame(name = input$file$name, datapath = input$file$datapath, stringsAsFactors = FALSE)
    })
    
    selectedSampleReactive <- reactive({ 1L })
    
    
    
    safeAnalyse <- function(path, guide, trim5, trim3){
      tryCatch(AnalyseSangerSample(readsangerseq(path), guide, trim5, trim3),
               error = function(e) list(error = conditionMessage(e)))
    }
    
    edited.analysesReactive <- reactive({
      files <- edited.filesReactive()
      guide <- guideReactive()
      trim5 <- input$trim5
      trim3 <- input$trim3
      lapply(files$datapath, safeAnalyse, guide = guide, trim5 = trim5, trim3 = trim3)
    })
    
    nt.Reactive <- reactive({
      req(!input$example, input$file_nt)
      safeAnalyse(input$file_nt$datapath, guideReactive(), input$trim5, input$trim3)
    })
    
    het_positionReactive <- reactive({
      if(!isTRUE(input$heterozygous) || input$example) return(NULL)
      pos <- input$het_position
      validate(need(!is.null(pos) && !is.na(pos) && pos >= 1 &&
                      pos <= length(guideReactive()) && pos == round(pos),
                    "Please enter a valid heterozygous position within the guide"))
      as.integer(pos)
    })
    
    expectedPositionReactive <- reactive({
      pos <- input$expected.edit.position
      if(is.null(pos) || length(pos) != 1 || is.na(pos)) return(NULL)
      if(pos < 1 || pos > length(guideReactive()) || pos != round(pos)) return(NULL)
      as.integer(pos)
    })
    
    sampleTableData <- function(a){
      EditingTableData(a$editing.df, a$sangs.filt, a$null.m.params, p.val.Reactive())
    }
    
    sampleFigure <- function(a, highlight, title){
      table.data <- sampleTableData(a)
      GuideFigure(GuideChromatogramPlot(GuideTraceData(a$basecalls, a$guide.coord),
                                        length(table.data$guide.seq), highlight),
                  GuideTablePlot(table.data, highlight), title)
    }
    
    normalisedFigure <- function(a, nt, het.position, title){
      norm.long <- NormaliseToWildType(a$editing.df, nt$editing.df, het.position)
      table.data <- NormalisedTableData(norm.long, sampleTableData(a), p.val.Reactive())
      GuideFigure(GuideChromatogramPlot(NormalisedTraceData(norm.long), length(table.data$guide.seq),
                                        het.position, y.label = "Normalised base area (%)", y.max = 100),
                  GuideTablePlot(table.data, het.position), title)
    }
    
    # PNG and SVG export of chromatogram + table as one figure
    figureDownloads <- function(id, label.fun, figure.reactive){
      output[[paste0("download_", id, "_png")]] <- downloadHandler(
        filename = function() paste0("EditR_", label.fun(), ".png"),
        content = function(file) SaveGuideFigure(figure.reactive(), file, "png"))
      output[[paste0("download_", id, "_svg")]] <- downloadHandler(
        filename = function() paste0("EditR_", label.fun(), ".svg"),
        content = function(file) SaveGuideFigure(figure.reactive(), file, "svg"))
    }
    
    sampleLabel <- function(name) tools::file_path_sans_ext(name)
    
    # outputs for edited sample i (edited and normalised figure + downloads), created once per index
    makeSampleOutputs <- function(i){
      local({
        edited.fig <- reactive({
          a <- edited.analysesReactive()[[i]]
          req(!is.null(a), is.null(a$error))
          sampleFigure(a, expectedPositionReactive(), edited.filesReactive()$name[i])
        })
        norm.fig <- reactive({
          a <- edited.analysesReactive()[[i]]
          nt <- nt.Reactive()
          het.position <- het_positionReactive()
          req(!is.null(a), is.null(a$error), is.null(nt$error), !is.null(het.position),
              nrow(a$editing.df) == nrow(nt$editing.df))
          normalisedFigure(a, nt, het.position,
                           paste0("Normalised base areas: ", edited.filesReactive()$name[i]))
        })
        qc.data <- reactive({
          a <- edited.analysesReactive()[[i]]
          req(!is.null(a), is.null(a$error))
          a
        })
        output[[paste0("qc_prefilter_", i)]] <- renderPlot({ qcPrefilterPlot(qc.data()) })
        output[[paste0("qc_signal_", i)]] <- renderPlotly({ qcSignalNoisePlot(qc.data()) })
        output[[paste0("qc_noise_", i)]] <- renderPlot({ qcNoisePercPlot(qc.data()) })
        output[[paste0("qc_chrom_", i)]] <- renderPlot({ qcChromatogram(qc.data()) })
        output[[paste0("baseinfo_table_", i)]] <- renderTable({ baseInfoTable(baseInfoFor(qc.data())) })
        output[[paste0("quad_plot_", i)]] <- renderPlot({ grid::grid.draw(quadPlotFor(qc.data())) })
        output[[paste0("figure_edited_", i)]] <- renderPlot({ edited.fig() })
        output[[paste0("figure_norm_", i)]] <- renderPlot({ norm.fig() })
        figureDownloads(paste0("edited_", i),
                        function() paste0(sampleLabel(edited.filesReactive()$name[i]), "_edited"), edited.fig)
        figureDownloads(paste0("norm_", i),
                        function() paste0(sampleLabel(edited.filesReactive()$name[i]), "_normalised"), norm.fig)
      })
    }
    
    outputs.created <- 0
    observe({
      req(input$example || !is.null(input$file))
      n <- nrow(edited.filesReactive())
      while(outputs.created < n){
        outputs.created <<- outputs.created + 1
        makeSampleOutputs(outputs.created)
      }
    })
    
    errorText <- function(msg) p(style = "color: #b94a48;", msg)
    figureBlock <- function(id, height = "480px"){
      tagList(plotOutput(paste0("figure_", id), height = height),
              downloadButton(paste0("download_", id, "_png"), "Download PNG"),
              downloadButton(paste0("download_", id, "_svg"), "Download SVG"),
              br(), br())
    }
    
    # QC / base info / quad plot blocks; id is the edited index or "nt" for the non-treated sample
    qcBlock <- function(id, name, a){
      if(!is.null(a$error)) return(tagList(h2(name), errorText(a$error), hr()))
      tagList(h2(name),
              plotOutput(paste0("qc_prefilter_", id)),
              plotlyOutput(paste0("qc_signal_", id)),
              plotOutput(paste0("qc_noise_", id)),
              plotOutput(paste0("qc_chrom_", id)),
              hr())
    }
    infoBlock <- function(id, name, a){
      head <- h3(paste0("Base info: ", name))
      if(!is.null(a$error)) tagList(head, errorText(a$error))
      else tagList(head, tableOutput(paste0("baseinfo_table_", id)))
    }
    quadBlock <- function(id, name, a){
      head <- h3(paste0("Quad plot: ", name))
      if(!is.null(a$error)) tagList(head, errorText(a$error))
      else tagList(head, plotOutput(paste0("quad_plot_", id)))
    }
    
    # edited samples followed by the non-treated sample (if uploaded)
    allSampleBlocks <- function(block, nt.label){
      files <- edited.filesReactive()
      res <- edited.analysesReactive()
      blocks <- lapply(seq_len(nrow(files)), function(i) block(i, files$name[i], res[[i]]))
      if(!isTRUE(input$example) && !is.null(input$file_nt)){
        blocks <- c(blocks, list(block("nt", paste0(input$file_nt$name, " (", nt.label, ")"), nt.Reactive())))
      }
      blocks
    }
    
    output$qc_sections <- renderUI({ allSampleBlocks(qcBlock, "non-treated") })
    output$baseinfo_sections <- renderUI({ allSampleBlocks(infoBlock, "non-treated") })
    output$quad_sections <- renderUI({ allSampleBlocks(quadBlock, "non-treated") })
    
    nt.qc.data <- reactive({
      nt <- nt.Reactive()
      req(is.null(nt$error))
      nt
    })
    output$qc_prefilter_nt <- renderPlot({ qcPrefilterPlot(nt.qc.data()) })
    output$qc_signal_nt <- renderPlotly({ qcSignalNoisePlot(nt.qc.data()) })
    output$qc_noise_nt <- renderPlot({ qcNoisePercPlot(nt.qc.data()) })
    output$qc_chrom_nt <- renderPlot({ qcChromatogram(nt.qc.data()) })
    output$baseinfo_table_nt <- renderTable({ baseInfoTable(baseInfoFor(nt.qc.data())) })
    output$quad_plot_nt <- renderPlot({ grid::grid.draw(quadPlotFor(nt.qc.data())) })
    
    # on-/off-target editing summary of the successfully analysed edited samples
    editingSummaryReactive <- reactive({
      pos <- expectedPositionReactive()
      if(is.null(pos)) return(NULL)
      files <- edited.filesReactive()
      res <- edited.analysesReactive()
      ok <- which(vapply(res, function(a) is.null(a$error), logical(1)))
      if(length(ok) == 0) return(NULL)
      
      # heterozygous sample: use the Abbott-normalised percentages instead of the raw ones
      het.position <- het_positionReactive()
      nt <- if(!is.null(het.position) && !is.null(input$file_nt)) nt.Reactive() else NULL
      use.norm <- !is.null(nt) && is.null(nt$error)
      if(use.norm) ok <- ok[vapply(res[ok], function(a) nrow(a$editing.df) == nrow(nt$editing.df), logical(1))]
      if(length(ok) == 0) return(NULL)
      
      perc.list <- lapply(res[ok], function(a){
        if(use.norm) NormalisedPercMatrix(a$editing.df, nt$editing.df, het.position) else RawPercMatrix(a$editing.df)
      })
      out <- EditingSummaryLines(perc.list, res[[ok[1]]]$editing.df$guide.seq, sampleLabel(files$name[ok]), pos)
      out$normalised <- use.norm
      out
    })
        output$editing_summary <- renderUI({
      if(is.null(expectedPositionReactive())){
        return(p(em("Enter the expected editing position in the guide to quantify on-target and off-target editing.")))
      }
      ed <- editingSummaryReactive()
      req(!is.null(ed))
      wellPanel(h4("Editing summary"), p(strong(ed$on.target)), p(ed$off.target),
                if(isTRUE(ed$normalised)) p(em("Values are normalised to the non-treated sample (Abbott formula at the heterozygous position).")))
    })

    output$edited_sections <- renderUI({
      files <- edited.filesReactive()
      res <- edited.analysesReactive()
      lapply(seq_len(nrow(files)), function(i){
        if(!is.null(res[[i]]$error)){
          tagList(h3(files$name[i]), errorText(res[[i]]$error))
        } else figureBlock(paste0("edited_", i))
      })
    })
    
    output$nt_section <- renderUI({
      if(is.null(input$file_nt)) return(p("No non-treated sample uploaded."))
      nt <- nt.Reactive()
      if(!is.null(nt$error)) return(tagList(h3(input$file_nt$name), errorText(nt$error)))
      figureBlock("nt")
    })
    
    nt.fig <- reactive({
      nt <- nt.Reactive()
      req(is.null(nt$error))
      sampleFigure(nt, NULL, input$file_nt$name)
    })
    output$figure_nt <- renderPlot({ nt.fig() })
    figureDownloads("nt", function() paste0(sampleLabel(input$file_nt$name), "_non-treated"), nt.fig)
    
    output$norm_sections <- renderUI({
      het.position <- het_positionReactive()
      files <- edited.filesReactive()
      res <- edited.analysesReactive()
      if(is.null(input$file_nt)){
        return(p("Please upload a non-treated sample to normalise to."))
      }
      nt <- nt.Reactive()
      lapply(seq_len(nrow(files)), function(i){
        a <- res[[i]]
        problem <- if(!is.null(nt$error)) paste0("Non-treated sample: ", nt$error)
        else if(!is.null(a$error)) a$error
        else if(nrow(a$editing.df) != nrow(nt$editing.df))
          "Guide regions of the non-treated and edited sample have different lengths"
        if(!is.null(problem)) tagList(h3(files$name[i]), errorText(problem))
        else figureBlock(paste0("norm_", i))
      })
    })
###########################################
    quadPlotFor <- function(a){
      editing.df <- a$editing.df
      null.m.params <- a$null.m.params
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
      
      arrangeGrob(p.a, p.c, p.g, p.t)
    }



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
    
    # PDF report with every sample (edited and non-treated). Drawn directly with the pdf device,
    # so no pandoc / LaTeX installation is needed.
    reportPage <- function(title, grob){
      grid::grid.newpage()
      grid::grid.draw(gridExtra::arrangeGrob(
        grob, top = grid::textGrob(title, x = 0.03, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 13)),
        padding = grid::unit(0.5, "cm")))
    }
    
    reportTextPage <- function(title, lines){
      grid::grid.newpage()
      grid::grid.text(title, x = 0.05, y = 0.95, just = c("left", "top"), gp = grid::gpar(fontface = "bold", fontsize = 16))
      wrapped <- unlist(lapply(lines, function(x) if(nchar(x) == 0) "" else strwrap(x, width = 95, exdent = 4)))
      grid::grid.text(paste(wrapped, collapse = "\n"), x = 0.05, y = 0.90, just = c("left", "top"),
                      gp = grid::gpar(fontsize = 10, lineheight = 1.3))
    }
    
    # figure occupies the top half of a fresh page
    printHalfPage <- function(p){
      grid::grid.newpage()
      print(p, vp = grid::viewport(x = 0.5, y = 1, width = 1, height = 0.5, just = c("center", "top")))
    }

    reportSamplePages <- function(name, role, a, het.position, nt = NULL){
      label <- paste0(name, " (", role, ")")
      if(!is.null(a$error)){
        reportTextPage(label, c("This sample could not be analysed:", a$error))
        return(invisible())
      }
      # QC
      qc <- arrangeGrob(qcPrefilterPlot(a), qcSignalNoiseGg(a), qcNoisePercPlot(a), ncol = 1)
      reportPage(paste0("Data QC: ", label), qc)
      # chromatogram + percent area table
      printHalfPage(sampleFigure(a, if(role == "edited") expectedPositionReactive() else NULL, name))
      # base info and quad plot
      info <- baseInfoTable(baseInfoFor(a))
      info[] <- lapply(info, function(x) if(is.numeric(x)) signif(x, 4) else x)
      tab <- gridExtra::tableGrob(info, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 9))
      reportPage(paste0("Base info and quad plot: ", label),
                 arrangeGrob(tab, quadPlotFor(a), ncol = 1, heights = c(1, 4)))
      # normalised
      if(!is.null(nt) && !is.null(het.position)){
        if(!is.null(nt$error)) {
          reportTextPage(paste0("Normalised: ", name), paste0("Non-treated sample: ", nt$error))
        } else if(nrow(a$editing.df) != nrow(nt$editing.df)) {
          reportTextPage(paste0("Normalised: ", name),
                         "Guide regions of the non-treated and edited sample have different lengths")
        } else {
          printHalfPage(normalisedFigure(a, nt, het.position, paste0("Normalised base areas: ", name)))
        }
      }
    }
    
    output$downloadReport <- downloadHandler(
      filename = function() paste0("EditR_report_", format(Sys.Date(), "%Y%m%d"), ".pdf"),
      content = function(file) {
        files <- edited.filesReactive()
        res <- edited.analysesReactive()
        het.position <- het_positionReactive()
        has.nt <- !isTRUE(input$example) && !is.null(input$file_nt)
        nt <- if(has.nt) nt.Reactive() else NULL
        
        trim <- if(!is.na(input$trim5) && !is.na(input$trim3)) paste0("Trimmed region: ", input$trim5, " to ", input$trim3)
        
        ed <- editingSummaryReactive()
        summary.lines <- if(is.null(ed)) character(0) else c("Editing summary (edited samples):", ed$on.target, ed$off.target,
                                                                        if(isTRUE(ed$normalised)) "Values are normalised to the non-treated sample (Abbott formula at the heterozygous position).")

        grDevices::pdf(file, width = 8.27, height = 11.69)
        on.exit(grDevices::dev.off(), add = TRUE)
        
        reportTextPage(paste0("EditR v", editrversion, " report"), c(
          paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
          paste0("Guide sequence: ", as.character(guideReactive())),
          paste0("P-value cutoff: ", p.val.Reactive()),
          paste0("Heterozygous position in guide: ", if(is.null(het.position)) "not used" else het.position),
          paste0("Expected editing position in guide: ", if(is.null(expectedPositionReactive())) "not set" else expectedPositionReactive()),
          trim,
          "",
          summary.lines,
          "",
          "Edited samples:",
          paste0("  - ", files$name),
          "",
          paste0("Non-treated sample: ", if(has.nt) input$file_nt$name else "none")
        ))
        
        for(i in seq_len(nrow(files))){
          reportSamplePages(files$name[i], "edited", res[[i]], het.position, nt)
        }
        if(has.nt) reportSamplePages(input$file_nt$name, "non-treated", nt, NULL)
      },
      contentType = "application/pdf"
    )

})
