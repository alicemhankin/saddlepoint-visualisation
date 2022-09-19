library(shiny)
library(shinydashboard)
library(mvtnorm)
library(numDeriv)
library(Deriv)
library(shinyWidgets)
library(RANN)
library(glue)
library(gstat)
library(plotly)
library(MASS)
library(scales)
options(shiny.useragg = TRUE)

scale_points = function(x1_points, x2_points, y_points, s, t, m) y_points*c(m)*exp(-s*x1_points-t*x2_points) 

extend_points_reactive = function(x1, x2, y, y_simulated, simulate_to_y, distr, vars) {
  if (simulate_to_y > y_simulated()) {
    new_count = rpois(1, simulate_to_y - y_simulated())
    new_x = distr()(new_count, vars)
    new_x1 = new_x[,1]
    new_x2 = new_x[,2]
    new_y = runif(new_count, min = y_simulated(), max = simulate_to_y)
    x1(c(x1(), new_x1))
    x2(c(x2(), new_x2))
    y(c(y(), new_y))
    y_simulated(simulate_to_y)
  }
}

ui <- dashboardPage(
  
  dashboardHeader(title = "A Visual Representation of the Bivariate Saddlepoint Approximation", titleWidth=700),
  
  
  # Sidebar to demonstrate various slider options ----
  dashboardSidebar(width=300,
                   
                   uiOutput("slidercols"),
                   
                   sliderInput("slider_val", "How many parameters do you need?",
                               min = 1, max = 20,
                               value = 2),
                   
                   textAreaInput("distrString", "R code for distribution as a function of integer n. 
                                 Requires a 2*n matrix output. Use vars[[1]] as your first parameter and vars[[2]] as your second, etc.",
                                 # value="function(n, vars){u1 = rgamma(n, vars[[1]])
                                 #        u2 = rgamma(n, shape=vars[[1]], rate=1)
                                 #        u12 = rgamma(n, vars[[2]])
                                 #        cbind(u1+u12, u2+u12)}",
                                 value="function(n, vars) rmvnorm(n, mean=c(vars[[1]], vars[[2]]),
                                 sigma=matrix(c(0.7,0.5,0.5,0.8), ncol=2))",
                                 width="100%", height="80px", resize="vertical"),
                   
                   textAreaInput("normalisingString", "R code for normalising factor as a function of s, t and vars:",
                                 # value="function(s,t,vars){exp(log(1-s)*(-vars[[1]]) + log(1-t)*(-vars[[1]]) + log(1-s-t)*(-vars[[2]]))}",
                                 value="function(s,t,vars) exp(t(c(s,t))%*%c(vars[[1]],vars[[2]]) + 
                                 0.5*t(c(s,t))%*%matrix(c(0.7,0.5,0.5,0.8), ncol=2)%*%c(s,t))",
                                 width="100%", height="100px", resize="vertical"),
                   
                   
                   fluidRow(
                     
                     column(6, sliderInput("slim", "s Range:",
                                           min = -10, max = 10,
                                           value = c(-5,5)),
                            
                            sliderInput("x1minmax", "x1 Range:",
                                        min=-20, max=20,
                                        value = c(-5,5))
                            , style='padding-right:0px;'),
                     column(6,
                            
                            sliderInput("tlim", "t Range:",
                                        min = -10, max = 10,
                                        value = c(-5,5)),
                            
                            
                            
                            sliderInput("x2minmax", "x2 Range:",
                                        min=-20, max=20,
                                        value = c(-5,5))
                            , style='padding-left:0px;')),
                   
                   actionButton("applydistr", "Apply", icon("paper-plane"),
                                style="color: #fff; background-color: var(--nav-color); border-color: #fff")
                   
                   
                   
  ),
  
  # Main panel for displaying outputs ----
  dashboardBody(
    
    #tags$head(tags$style(HTML('
    #  .main-header .sidebar-toggle:before {
    #    content: \' \\43 \\2009 \\4c \\2009 \\49 \\2009 \\43 \\2009 \\4b \\2009 \\21  \';}'))),
    uiOutput("changesidebartoggleicon"),
    
    #changing colour of active tab indicator reactively
    tags$style(".nav-tabs-custom .nav-tabs li.active {border-top-color: var(--nav-color);}"),
    
    #changing colour of header reactively
    uiOutput("navcolor"),
    tags$head(tags$style(HTML('
      .skin-blue .main-header .logo {background-color: var(--nav-color);}
      .skin-blue .main-header .logo:hover {background-color: var(--nav-color);}
      .skin-blue .main-header .navbar .sidebar-toggle:hover{background-color: var(--nav-color);}
      .skin-blue .main-header .navbar {background-color: var(--nav-color);}'))),
    
    
    # Output: plot summarizing the values entered ----
    fluidRow(box(width=3, 
                 
                 uiOutput("dynamic_widget"),
                 
                 hr(style = "border: 1px solid;"),
                 
                 radioGroupButtons(
                   inputId = "colourswitch",
                   choices = c( "Choose s,t"="false", "Choose K'(s,t)"="true"),
                   justified=TRUE
                 ),
                 
                 fluidRow(column(6, uiOutput("sLoopSlider"),
                                 uiOutput("x1LoopSlider")),
                          
                          column(6,uiOutput("tLoopSlider"),
                                 uiOutput("x2LoopSlider"))),
                 
                 prettyCheckbox("togglenormalising", "Normalise", TRUE)),
             
             tabBox(id="tabs", width=6, height="83vh", tabPanel("Saddlepoint Plot", plotOutput("scaled_plot", height="80vh")),
                    tabPanel("Tilted Density Plot", plotOutput("plot2", height="80vh")),
                    #tabPanel("3", plotlyOutput("plot3"))
                    ),
             
             box(width=3, 
                 
                 uiOutput("text"),
                 
                 uiOutput("density"),
                 
                 sliderInput("alpha1", "Heatmap Opacity:",
                             min = 0, max = 1, value = 0.7),
                 sliderInput("alpha2", "Point Opacity:",
                             min = 0, max = 1, value = 1),
                 
                 
                 uiOutput("x1slider"),
                 uiOutput("x2slider"),
                 
                 # Input: Specification of range within an interval ----
                 sliderInput("ylim", "y Range:",
                             min = 0, max = 1000,
                             value = c(0,100)),
                 sliderInput("IntensityToSimulate", "Intensity to simulate:",
                             min = 0, max = 100000, value = 100),
                 
                 prettyCheckbox("togglex", "Toggle Crosshairs", TRUE)
             ))
    
    # tags$div(style="margin-bottom:400px;", fluidRow(
    #tabBox(width=9,
    #tabPanel("Saddlepoint Approximation", plotOutput("plot1")),
    # tabPanel("Tilted Distribution", plotOutput("plot2")),
    #)
    # ))
  )
)


# Define server logic ----
server <- function(input, output) {
  
  
  
  x1 <- reactiveVal()  
  x2 <- reactiveVal()
  y <- reactiveVal()
  
  y_simulated <- reactiveVal(0)
  
  vars = reactive({
    a = c(rep(list(1),20))
    a[1:as.numeric(input$slider_val)] = lapply(1:as.numeric(input$slider_val),
                                               function(i) ifelse(is.null(input[[paste0("var",i)]]),1,input[[paste0("var",i)]]))
    a
  })
  
  
  distr <- reactiveVal(eval(parse(text=isolate(input$distrString))))
  observeEvent(c(input$applydistr,vars()),{
    distr(eval(parse(text=isolate(input$distrString))))
    x1(c());x2(c())
    y(c())
    y_simulated(0)
    extend_points_reactive(x1, x2, y, y_simulated, 
                           input$IntensityToSimulate, 
                           distr, vars())
  })
  normalising <- reactiveVal(eval(parse(text=isolate(input$normalisingString))))
  
  observeEvent(input$applydistr, {
    normalising(eval(parse(text=isolate(input$normalisingString))))
  })
  
  observeEvent(c(input$IntensityToSimulate,vars()), {
    extend_points_reactive(x1, x2, y, y_simulated, 
                           input$IntensityToSimulate, 
                           distr, vars())
  })
  
  m = reactive({if (input$togglenormalising) {
    normalising()(input$s, input$t, vars())
  } else {
    1
  }})
  
  scaled_y <- reactive({
    if(input$colourswitch=="false"){
      scale_points(x1(),x2(), y(), input$s, input$t, m())
    } else{
      st = get_k_prime_inverse(input$x1,input$x2)
      scale_points(x1(),x2(), y(), st[1], st[2], m())
    }
  })
  
  
  cgfs = reactive({
    m = normalising()
    
    ms = Deriv(m,"s")
    mt = Deriv(m, "t")
    
    
    mss = Deriv(ms, "s")
    mst = Deriv(ms, "t")
    mtt = Deriv(mt, "t")
    
    list("m" = m, "ms" = ms, "mt"=mt, "mss" = mss, "mst" = mst, "mtt" = mtt)  
  })
  
  get_ks = function(s,t,cc,vars){
    m = cc$m(s,t,vars)
    mt = cc$mt(s,t,vars)
    ms = cc$ms(s,t,vars)
    
    k = log(m)
    
    k_prime = c(ms / m, mt / m)
    
    k_dblprime_ss =  (cc$mss(s,t,vars)*m - ms^2) / m^2
    k_dblprime_tt =  (cc$mtt(s,t,vars)*m - mt^2) / m^2
    k_dblprime_st = (cc$mst(s,t,vars)*m - mt*ms) / m^2
    k_dblprime = matrix(c(k_dblprime_ss,k_dblprime_st,k_dblprime_st,k_dblprime_tt), byrow=T, nrow=2, ncol=2)
    
    list(k,k_prime,k_dblprime)
  }
  
  get_approx = function(s,t,cc,vars){
    ks = get_ks(s,t,cc,vars)
    exp(ks[[1]] - c(s,t)%*%ks[[2]]) /sqrt(det(ks[[3]] * 2 * pi))
  }
  
  xyz <-reactive({
    req(input$normalisingString)
    
    cc = cgfs()
    vars = vars()
    x1seq = seq(minmax_x()[1],minmax_x()[2],0.2)
    x2seq = seq(minmax_x()[3],minmax_x()[4],0.2)
    sseq = seq(st_lims()[1],st_lims()[2],0.2)
    tseq = seq(st_lims()[3],st_lims()[4],0.2)
    
    corresponding_t = rep(tseq, length.out = length(sseq)*length(tseq))
    corresponding_s = rep(sseq, each=length(tseq))
    
    x_for_each_st = mapply(function(s, t) get_ks(s, t, cc, vars)[[2]], corresponding_s, corresponding_t)
    x_for_each_st[x_for_each_st == Inf | x_for_each_st == -Inf | is.na(x_for_each_st)] = -9999999
    x_for_each_st = t(x_for_each_st)
    
    each_x = cbind(rep(x1seq, 1, each=length(x2seq)), rep(x2seq, length(x1seq)))
    indexes = nn2(x_for_each_st,each_x, 1)$nn.idx 
    
    s = corresponding_s[indexes]
    t = corresponding_t[indexes]
    
    approx = mapply(function(s, t) get_approx(s, t, cc, vars), s, t)
    
    cbind(each_x, approx)
    
  })  
  
  get_k_prime_inverse = function(x1, x2){
    req(input$normalisingString)
    
    cc = cgfs()
    vars = vars()
    sseq = seq(st_lims()[1],st_lims()[2],0.2)
    tseq = seq(st_lims()[3],st_lims()[4],0.2)
    
    corresponding_t = rep(tseq, length.out = length(sseq)*length(tseq))
    corresponding_s = rep(sseq, each=length(tseq))
    
    x_for_each_st = mapply(function(s, t) get_ks(s, t, cc, vars)[[2]], corresponding_s, corresponding_t)
    x_for_each_st[x_for_each_st == Inf | x_for_each_st == -Inf | is.na(x_for_each_st)] = -9999999
    x_for_each_st = t(x_for_each_st)
    
    x = cbind(x1, x2)
    index = nn2(x_for_each_st,x, 1)$nn.idx 
    
    s = corresponding_s[index]
    t = corresponding_t[index]
    
    c(s,t)
  }
  
  
  x1x2 = reactive({
    if(input$colourswitch=="false"){
      get_ks(input$s, input$t, cgfs() ,vars())
    }else{
      st = get_k_prime_inverse(input$x1,input$x2)
      get_ks(st[1], st[2], cgfs() ,vars())
    }
  })
  
  
  output$scaled_plot <- renderPlot({
    req(input$s, input$t, y_simulated())
    
    title(paste("s = ", input$s, "t=", input$t))
    
    par(mar=rep(0,4))
    layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE),
           widths=c(1, lcm(5)),heights=c(1, lcm(5)))
    cols = c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")[c(2,3,4,6,7)]
    
    xyz = xyz()
    one  = xyz2img(xyz, tolerance=100*.Machine$double.eps)
    two = cbind(rep(1,length(unique(xyz[,2]))), aggregate(xyz, list(xyz[,2]), FUN=sum)[,c(1,4)])
    three = cbind(rep(1,length(unique(xyz[,1]))), aggregate(xyz, list(xyz[,1]), FUN=sum)[,c(1,4)])
    
    indices = (y() <= min(input$ylim[2],input$IntensityToSimulate))
    
    x1 = ifelse(is.na(x1()[indices]), 9999, x1()[indices])
    x2 = ifelse(is.na(x2()[indices]), 9999, x2()[indices])
    x1x2=x1x2()
    
    
    pointdensity = kde2d(x1()[indices],x2()[indices], 
                         n = c(length(seq(input$x1lim[1],input$x1lim[2],0.2)), length(seq(input$x2lim[1],input$x2lim[2],0.2))), 
                         lims = c(input$x1lim, input$x2lim))
    
    if(input$densityswitch){
      image(one,
            col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1), 
            xlim = c(input$x1lim[1],input$x1lim[2]), ylim =c(input$x2lim[1],input$x2lim[2]))
    }else{
      image(pointdensity, 
            col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1), 
            xlim = c(input$x1lim[1],input$x1lim[2]), ylim =c(input$x2lim[1],input$x2lim[2]))
    }
    
    points(x1()[indices],x2()[indices], xlim=input$x1lim, ylim=input$x2lim, 
           col=scales::alpha("black",input$alpha2),
           pch=16)
    
    if(input$togglex){
    abline(h=x1x2[[2]][2])
    abline(v=x1x2[[2]][1])
    }
    
    
    if(input$densityswitch){
      image(y=two[,2], x=1, z=matrix(two[,3], nrow=1),
            col= hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1),ylim=input$x2lim, xlim=input$ylim)
    }else{
      image(rbind(colSums(pointdensity[[3]])),ylim=input$x2lim,y=seq(input$x2lim[1],input$x2lim[2],0.2),xlim=input$ylim,
            col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1))
    }
    
    points(scaled_y()[indices], x2()[indices], ylim=input$x2lim, xlim=input$ylim, 
           col=scales::alpha("black",input$alpha2),
           pch=16)
    
    if(input$togglex){abline(h=x1x2[[2]][2])}
    
    
    if(input$densityswitch){
      image(x=three[,2], y=1, z=matrix(three[,3], ncol=1), xlim=input$x1lim, ylim=input$ylim,
            col= hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1))
    }else{
      image(z=cbind(rowSums(pointdensity[[3]])),x=seq(input$x1lim[1],input$x1lim[2],0.2),ylim=input$ylim,xlim=input$x1lim,
            col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1))
    }
    points(x1()[indices], scaled_y()[indices], xlim=input$x1lim, ylim=input$ylim,
           col=scales::alpha("black",input$alpha2),
           pch=16)
    if(input$togglex){abline(v=x1x2[[2]][1])}
    
    
    
    
  })
  
  output$plot2 <- renderPlot({
    req(input$normalisingString)
    req(input$s, input$t, y_simulated(),input$x1,input$x2)
    
    par(mar=rep(0,4))
    layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE),
           widths=c(1, lcm(5)),heights=c(1, lcm(5)))
    
    x1seq = seq(minmax_x()[1],minmax_x()[2],0.2)
    x2seq = seq(minmax_x()[3],minmax_x()[4],0.2)
    indices = (scaled_y() <= input$ylim[2])
    x1x2=x1x2()
    
    mat = matrix(rep(NA, length(x1seq)*length(x2seq)), nrow=length(x1seq))
    
    for (i in 1:length(x1seq)){
      for (j in 1:length(x2seq)){
        mat[i,j] = Rfast::dmvnorm(x=c(x1seq[i],x2seq[j]), mu = x1x2[[2]], 
                                  sigma =x1x2[[3]])
      }
    }
    
    pointdensity = kde2d(x1()[indices],x2()[indices], 
                         n = c(length(seq(input$x1lim[1],input$x1lim[2],0.2)), length(seq(input$x2lim[1],input$x2lim[2],0.2))), 
                         lims = c(input$x1lim, input$x2lim))
    
    
    #zlims = c(0,max(pointdensity[[3]], mat))


    
    
    if(input$densityswitch){
      image(mat, x=x1seq, y=x2seq, xlim = c(input$x1lim[1],input$x1lim[2]), ylim =c(input$x2lim[1],input$x2lim[2]),
          col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1))
    }else{
      image(pointdensity, 
            col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1), 
            xlim = c(input$x1lim[1],input$x1lim[2]), ylim =c(input$x2lim[1],input$x2lim[2]))
    }
    points(x1()[indices], x2()[indices], xlim=input$x1lim, ylim=input$x2lim, 
           col=scales::alpha("black",input$alpha2),
           pch=16)
    if(input$togglex){
    abline(h=x1x2[[2]][2])
    abline(v=x1x2[[2]][1])
    }
    
    
    if(input$densityswitch){
      image(rbind(colSums(mat)),ylim=input$x2lim,y=x2seq,xlim=input$ylim,
          col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1))
    }else{
      image(rbind(colSums(pointdensity[[3]])),ylim=input$x2lim,y=seq(input$x2lim[1],input$x2lim[2],0.2),xlim=input$ylim,
            col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1))
    }
    
    
    points(scaled_y()[indices], x2()[indices], ylim=input$x2lim, xlim=input$ylim, col=scales::alpha("black",input$alpha2),
           pch=16)
    if(input$togglex){abline(h=x1x2[[2]][2])}
    
    if(input$densityswitch){
      image(z=cbind(rowSums(mat)),x=x1seq,ylim=input$ylim,xlim=input$x1lim,
          col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1))
    }else{
      image(z=cbind(rowSums(pointdensity[[3]])),x=seq(input$x1lim[1],input$x1lim[2],0.2),ylim=input$ylim,xlim=input$x1lim,
            col = hcl.colors(200, "Blues3", rev = TRUE, alpha=input$alpha1))
    }
    
    
    points(x1()[indices], scaled_y()[indices], xlim=input$x1lim, ylim=input$ylim, col=scales::alpha("black",input$alpha2),
           pch=16)
    if(input$togglex){abline(v=x1x2[[2]][1])}
    
    # a = c()
    # for(i in x1seq){
    #   x2 = max(input$IntensityToSimulate*c(normalising()(input$s, input$t, vars()))
    #                  * exp(-input$s*i -input$t*x2()[indices]))
    #   a = c(a, x2)
    # }
    # lines(x1seq,a, col="red")
    
  })
  
  output$plot3 <- renderPlotly({
    x1seq = seq(minmax_x()[1],minmax_x()[2],0.2)
    x2seq = seq(minmax_x()[3],minmax_x()[4],0.2)
    indices = (scaled_y() <= input$ylim[2])
    
    df = data.frame(cbind(x1()[indices], x2()[indices], scaled_y()[indices]))
    colnames(df)=c("x1", "x2", "y")
    plot_ly(df, x = ~x1, y = ~x2, z=~y, type="scatter3d")
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  st_lims <- reactiveVal(c(isolate(input$slim[1]),isolate(input$slim[2]),isolate(input$tlim[1]),isolate(input$tlim[2])))
  
  observeEvent(input$applydistr, {st_lims(c(input$slim[1],input$slim[2],input$tlim[1],input$tlim[2]))})
  
  
  output$sLoopSlider <- renderUI({
    conditionalPanel(condition = "input.colourswitch == 'false'",
                     sliderInput("s", "s:",
                                 min=st_lims()[1], max=st_lims()[2],
                                 value = 0, step = 0.1,
                                 animate =
                                   animationOptions(interval = 300, loop = TRUE)))
  })
  
  output$tLoopSlider <- renderUI({
    conditionalPanel(condition = "input.colourswitch == 'false'",
                     sliderInput("t", "t:",
                                 min=st_lims()[3], max=st_lims()[4],
                                 value = 0, step = 0.1,
                                 animate = animationOptions(interval = 300, loop = TRUE)))
  })
  
  output$x1LoopSlider <- renderUI({
    conditionalPanel(condition = "input.colourswitch == 'true'",
                     sliderInput("x1", "Desired mean x1 value:",
                                 min=minmax_x()[1], max=minmax_x()[2],
                                 value = 0, step = 0.1,
                                 animate = animationOptions(interval = 300, loop = TRUE)))
  })
  
  output$x2LoopSlider <- renderUI({
    conditionalPanel(condition = "input.colourswitch == 'true'",
                     sliderInput("x2", "Desired mean x2 value:",
                                 min=minmax_x()[3], max=minmax_x()[4],
                                 value = 0, step = 0.1,
                                 animate = animationOptions(interval = 300, loop = TRUE)))
  })
  
  output$dynamic_widget <- renderUI({
    num <- as.integer(input$slider_val)
    column(12,lapply(1:num,function(i) {
      column(4,numericInput(inputId = paste0("var",i),label=paste0("vars[[",i,"]]"),
                            value=0.5,step=0.1), style='padding:2px;')
    }), style='padding:0px;')
  })
  
  minmax_x <- reactiveVal(c(isolate(input$x1minmax)[1],isolate(input$x1minmax)[2], isolate(input$x2minmax)[1],
                            isolate(input$x2minmax)[2]))
  
  observeEvent(input$applydistr, {minmax_x(c(input$x1minmax[1],input$x1minmax[2], input$x2minmax[1], input$x2minmax[2]))})
  
  output$x1slider = renderUI({sliderInput("x1lim", "x1 zoom:",
                                          min = minmax_x()[1], max = minmax_x()[2],
                                          value = c(-3,3))})
  
  output$x2slider = renderUI({sliderInput("x2lim", "x2 zoom:",
                                          min = minmax_x()[3], max = minmax_x()[4],
                                          value = c(-3,3))})
  
  output$slidercols <- renderUI({
    colour <- "#CC79A7"
    tagList(
      tags$style(sprintf(".irs-bar-edge, .irs-bar, .irs-single, .irs-from, .irs-to {background: %s !important;}", colour)),
      tags$style(sprintf(".irs-bar {border-top: %s !important;}", colour)),
      tags$style(sprintf(".irs-bar {border-bottom: %s !important;}", colour))
    )
  })
  
  n <- 0
  makeReactiveBinding('n')
  observeEvent(list(input$sidebarExpanded,input$sidebarCollapsed), { 
    n <<- n + 1
  })
  
  output$changesidebartoggleicon <- renderUI({
    if(n %% 2 == 1){
      tags$head(tags$style(HTML('
      .main-header .sidebar-toggle:before {
        content: \' \\f060 \';}')))
    }else{
      tags$head(tags$style(HTML('
      .main-header .sidebar-toggle:before {
        content: \' \\f061 \';}')))
    }
  })
  
  
  output$navcolor <- renderUI({
    tags$style(HTML(glue(
      ":root {--nav-color: @{'#CC79A7'}@}",
      .open = "@{", .close = "}@"
    )))
  })
  
  output$text <- renderUI({HTML("<b>What should the blue heatmap show?</b>")})
  
  output$density <- renderUI({
    if(input$tabs == "Saddlepoint Plot"){
      radioGroupButtons(
        inputId = "densityswitch",
        choiceValues = c(T,F),
        choiceNames = c(HTML("Saddlepoint<br>Approximation"), HTML("Point<br>Density")),
        justified=TRUE)
    }else{
      radioGroupButtons(
        inputId = "densityswitch",
        choiceValues = c(T,F),
        choiceNames = c(HTML("Normal Approx. to<br>the Tilted Density"), HTML("Point<br>Density")),
        justified=TRUE)
    }
  })
  

  
}

# Create Shiny app ----
shinyApp(ui, server)
