library(shinydashboard)
library(shiny)
library(numDeriv)
library(latex2exp)
library(shinyWidgets)
library(fresh)
library(bslib)
library(shinyjs)
library(shinyBS)
library(glue)
library(Deriv)
library(plyr)
options(shiny.useragg = TRUE)     #https://github.com/rstudio/shiny/issues/3626   #https://github.com/s-u/Cairo/issues/37
#setwd("C:/Users/alice/Desktop/793/")

distfunc = function(distrString){
  
  names = c("ddist", "rdist", "m", "firstminmax", "secondminmax", "step", "max_s", "slidernames", "scalex", "discrete",
            "condition", "min_x")
  mylist <- sapply(names,function(x) NULL)
  
  def_step = 0.05
  mylist$min_x = -10
  
  switch(distrString,
         
         "Normal" = {mylist$start = c(0,1);
         mylist$ddist = function(x, var1, var2) dnorm(x, var1, var2);
         mylist$rdist = function(n) rnorm(n, mean=mylist$start[1],
                                          sd = mylist$start[2]);
         mylist$m = function(s, var1, var2) exp(s*var1+0.5*s^2*var2^2);
         mylist$firstminmax = c(-10,10);
         mylist$secondminmax = function(var1, var2) c(def_step,3.4);
         mylist$step = rep(def_step, 2);
         mylist$max_s = function(var1, var2) 10;
         mylist$slidernames = c("Mean:", "Standard Deviation:");
         mylist$scalex = function(new1,new2,x) {a = new2/mylist$start[2];
         b = new1 - a * mylist$start[1];
         a * x + b};
         mylist$discrete = F
         mylist$start = c(0,1)
         },
         
         "Exponential"= {mylist$start = c(2,0);
         mylist$ddist = function(x, var1, var2) ifelse(x<0,NA,dexp(x, var1));
         mylist$rdist = function(n) rexp(n, rate=mylist$start[1]);
         mylist$m = function(s, var1, var2) 1/(1-s*(1/var1));
         mylist$firstminmax = c(def_step,10);
         mylist$secondminmax = function(var1, var2) c(-10,10);
         mylist$step = rep(def_step, 2);
         mylist$max_s = function(var1, var2) max(0,var1 - def_step);
         mylist$slidernames = c("Rate:", "");
         mylist$scalex = function(new1,new2,x) mylist$start[1] * x / new1;
         mylist$discrete = F;
         mylist$min_x = def_step
         },
         
         "Geometric" = {mylist$start = c(0.5,0);
         mylist$ddist = function(x, var1, var2) (1-var1)^x*var1;
         mylist$rdist = function(n) runif(n);
         mylist$m = function(s, var1, var2) var1/(1-exp(s)+var1*exp(s));
         mylist$firstminmax = c(def_step,1-def_step);
         mylist$secondminmax = function(var1, var2) c(-10,10);
         mylist$step = rep(def_step,2);
         mylist$max_s = function(var1, var2) max(0,floor(-suppressWarnings(log(1-var1))*100)/100, na.rm=T);
         mylist$slidernames = c("Probability of Success:", "");
         mylist$scalex = function(new1,new2,x) suppressWarnings(qgeom(x, new1));
         mylist$discrete = T
         mylist$min_x = def_step
         },
         
         "Gamma" = {mylist$start = c(1,1);
         mylist$ddist = function(x, var1, var2) ifelse(x<0,NA,dgamma(x, var1, var2));
         mylist$rdist = function(n) rgamma(n,shape=mylist$start[1],rate=mylist$start[2]);
         mylist$m = function(s, var1, var2) (1-(s/var2))^(-1*var1);
         mylist$firstminmax = c(def_step,10);
         mylist$secondminmax = function(var1, var2) c(def_step,10);
         mylist$step = rep(def_step, 2);
         mylist$max_s = function(var1, var2) max(0,var2 - def_step);
         mylist$slidernames = c("Shape:", "Rate:");
         mylist$scalex = function(new1,new2,x) suppressWarnings(qgamma(pgamma(x, mylist$start[1], mylist$start[2]), new1, new2));
         mylist$discrete = F
         mylist$min_x = def_step
         },
         
         "Chi-Squared" = {mylist$start = c(2,0);
         mylist$ddist = function(x, var1, var2) ifelse(x<0,NA,dchisq(x, var1));
         mylist$rdist = function(n) rchisq(n,df=mylist$start[1]);
         mylist$m = function(s, var1, var2) (1-2*s)^(-0.5*var1);
         mylist$firstminmax = c(def_step,10);
         mylist$secondminmax = function(var1, var2) c(-10,10);
         mylist$step = rep(def_step,2);
         mylist$max_s = function(var1, var2) 0.5-def_step;
         mylist$slidernames = c("Degrees of Freedom:", "");
         mylist$scalex = function(new1,new2,x) qchisq(pchisq(x,mylist$start[1]),new1);
         mylist$discrete = F
         mylist$min_x = def_step
         },
         
         "Poisson" = {mylist$start = c(2,0);
         mylist$ddist = function(x, var1, var2) ifelse(x < -1,NA,suppressWarnings(var1^x * exp(-var1) / factorial(x)));
         mylist$rdist = function(n) runif(n);
         mylist$m = function(s, var1, var2) exp(var1*(expm1(s)));
         mylist$firstminmax = c(def_step,10);
         mylist$secondminmax = function(var1, var2) c(-10,10);
         mylist$step = rep(def_step, 2);
         mylist$max_s = function(var1, var2) 6.2 - log(var1, 3); #this was essentially found through trial and error. not sure what it should be
         mylist$slidernames = c("Rate:", "");
         mylist$scalex = function(new1,new2,x) suppressWarnings(qpois(x, new1));
         mylist$discrete = T;
         mylist$min_x = 0.15
         },
         
         "Binomial" = {mylist$start = c(0.5,10);
         mylist$ddist = function(x, var1, var2) ifelse(x>var2+1,NA,suppressWarnings(factorial(var2) * var1^x * (1-var1)^(var2-x)/(factorial(x)*factorial(var2-x))));
         mylist$rdist = function(n) runif(n);
         mylist$m = function(s, var1, var2) (1-var1+var1*exp(s))^var2;
         mylist$firstminmax = c(def_step,1-def_step);
         mylist$secondminmax = function(var1, var2) c(1,10);
         mylist$step = c(def_step,1);
         mylist$max_s = function(var1, var2) var2;
         mylist$slidernames = c("Probability of a success:", "Number of Trials:");
         mylist$scalex = function(new1,new2,x) suppressWarnings(qbinom(x, new2, new1));
         mylist$discrete = T
         mylist$min_x = -1
         },
         
         "Negative Binomial" = {mylist$start = c(0.5,2);
         mylist$ddist = function(x, var1, var2) ifelse(x < -1,NA,suppressWarnings(factorial(x+var2-1) * var1^var2 * (1-var1)^x / (factorial(x) * factorial(var2-1))));
         mylist$rdist = function(n) runif(n);
         mylist$m = function(s, var1, var2) (var1/(1-exp(s)+var1*exp(s)))^var2;
         mylist$firstminmax = c(def_step,1);
         mylist$secondminmax = function(var1, var2) c(1,10);
         mylist$step = c(def_step,1);
         mylist$max_s = function(var1, var2) -suppressWarnings(log(1-var1));
         mylist$slidernames = c("Probability of success:", "Target number of successes:");
         mylist$scalex = function(new1,new2,x) suppressWarnings(qnbinom(x, new2, new1));
         mylist$discrete = T
         mylist$min_x = 0.3
         },
         
         "Uniform" = {mylist$start = c(0,2);
         mylist$ddist = function(x, var1, var2) suppressWarnings(dunif(x, var1, var2));
         mylist$rdist = function(n) runif(n, mylist$start[1], mylist$start[2]);
         mylist$m = function(s, var1, var2) ifelse(s == 0, 1, (exp(var1*s) * expm1((var2-var1)*s)) / (s * (var2-var1)));
         mylist$firstminmax = c(-10,10);
         mylist$secondminmax = function(var1, var2) c(var1+def_step,10);
         mylist$step = rep(def_step, 2);
         mylist$max_s = function(var1, var2) 10;
         mylist$slidernames = c("Minimum:", "Maximum:");
         mylist$scalex = function(new1,new2,x) suppressWarnings(qunif(punif(x,mylist$start[1],mylist$start[2]),new1,new2));
         mylist$discrete = F
         }
         
  )
  
  mylist$condition = ifelse(mylist$slidernames[2]=="", "false", "true")
  
  mylist
}

scale_points_y = function(x_points, y_points, s, m) y_points*m*exp(-s*x_points)

extend_points_reactive = function(x1, x2, y, y_simulated, simulate_to_y, distr1, distr2) {
  if (simulate_to_y > y_simulated()) {
    new_count = rpois(1, simulate_to_y - y_simulated())
    new_x1 = distr1()(new_count)
    new_x2 = distr2()(new_count)
    new_y = runif(new_count, min = y_simulated(), max = simulate_to_y)
    x1(c(x1(), new_x1))
    x2(c(x2(), new_x2))
    y(c(y(), new_y))
    y_simulated(simulate_to_y)
  }
}


ui <- dashboardPage(
  
  
  dashboardHeader(
    title = "A Visual Representation of the Saddlepoint Approximation",
    titleWidth = 600),
  
  
  dashboardSidebar(
    
    width=420,
    
    fluidRow(
      column(6, selectInput("distrString1", label = "Distribution for x1:", 
                            choices = c("Normal", "Exponential","Gamma", "Chi-Squared", "Poisson",
                                        "Geometric", "Negative Binomial", "Binomial"), width="100%"),
             
             htmlOutput("name_first_slider"),
             
             uiOutput("first_slider"),
             
             htmlOutput("name_second_slider"),
             
             uiOutput("second_slider"),
             
             style = "padding-right:0px; margin-right:-5px"),
      
      column(6,
             
             selectInput("distrString2", label = "Distribution for x2:", 
                         choices = c("Normal", "Exponential","Gamma", "Chi-Squared", "Poisson",
                                     "Geometric", "Negative Binomial", "Binomial"), width="100%"),
             
             htmlOutput("name_third_slider"),
             
             uiOutput("third_slider"),
             
             htmlOutput("name_fourth_slider"),
             
             uiOutput("fourth_slider"),
             
             style = "padding-left:0px; margin-left:-5px;padding-right:5px; margin-right:0px;")),
    
    actionButton("applydistr", "Re-Simulate", icon("refresh", lib="glyphicon"),
                 style="color: #FFFFFF; background-color: var(--nav-color)"),
    
    hr(),
    
    column(12, radioGroupButtons(
      inputId = "colourswitch",
      choices = c( "Choose s"="false", "Choose K'(s)"="true"),
      justified=TRUE
    ),
    
    uiOutput("xLoopSlider"),
    
    
    uiOutput("sLoopSlider"),
    
    fluidRow(column(6,
                    
                    uiOutput("x_limit"),
                    
                    prettyCheckbox("togglecol", "Toggle Colourblind Mode", FALSE),
                    
                    style = "padding-right:0px; margin-right:-5px; margin-top:-10px;"
                    
    ),
    
    column(6,sliderInput("IntensityToSimulate", "Intensity to simulate:",
                         min = 0, max = 100000, value = 10000),
           style = "padding-left:0px; padding-right:5px; margin-left:-5px; margin-top:-10px;"
    )),
    
    style = "margin-top:-20px; padding:0px")
    
    
    # HTML("<hr style='height:50px; visibility:hidden; margin-bottom:-1px;' />"),
  ),
  
  
  dashboardBody(
    
    tags$head(tags$style(HTML('
      .main-header .sidebar-toggle:before {
        content: \' \\43 \\2009 \\4c \\2009 \\49 \\2009 \\43 \\2009 \\4b \\2009 \\21  \';}'))),
    
    
    #Changing colour of radio button labels (toggle graph lines on/off)
    tags$style(".tilt {color: #E69F00;}"),
    tags$style(".norm {color: #0072B2;}"),
    tags$style(".sim {color: #000000;}"),
    tags$style(".sadl {color: #0072B2;}"),
    tags$style(".true {color: #E69F00;}"),
    
    #defines nav-color
    tags$style(HTML(glue(
      ":root {--nav-color: @{'#E69F00'}@}",
      .open = "@{", .close = "}@"
    ))),
    
    #changing colour of active tab indicator
    tags$style(".nav-tabs-custom .nav-tabs li.active {border-top-color: var(--nav-color);}"),
    
    #changing colour of navbar
    tags$head(tags$style(HTML('
      .skin-blue .main-header .logo {background-color: var(--nav-color);}
      .skin-blue .main-header .logo:hover {background-color: var(--nav-color);}
      .skin-blue .main-header .navbar .sidebar-toggle:hover{background-color: var(--nav-color);}
      .skin-blue .main-header .navbar {background-color: var(--nav-color);}'))),
    
    
    uiOutput("changesidebartoggleicon"),
    
    uiOutput("slidercols"),
    
    #setting background colour
    setBackgroundColor("#efefef", shinydashboard = T),
    
    fluidRow(
      
      box(width=9,plotOutput("scaled_plot_old", height = "38vh"), height = "41vh"),
      
      box(width=3,title = "Controls",
          sliderInput("ylim", "y Range:",min = 0, max = 1000,value = c(0,100)),
          prettyCheckbox("togglex2", "Toggle Crosshairs", TRUE),
          prettyCheckbox("togglenormalising", "Normalise", TRUE))
    ),
    
    tags$div(style="margin-bottom:1px;",
             
             fluidRow(
               tabBox(width=9,id="tabs",
                      tabPanel("Saddlepoint Approximation", value="sadd", plotOutput("another_plot_old", height = "38vh")),
                      tabPanel("Tilted Distribution",  value="tilt",plotOutput("yet_another_plot_old", height = "38vh"))
                      , height = "41vh"),
               
               box(width=3, title = "Controls",
                   sliderInput("ylim2", "y Range:", min = 0, max = 2, value = c(0,0.5), step=0.01),
                   prettyCheckbox("togglelog", "Toggle Log", FALSE),
                   prettyCheckbox("togglex", "Toggle Crosshairs", TRUE)
               ),
               
               conditionalPanel(
                 condition = "input.tabs == 'tilt'",
                 box(width=3, tags$div(class="norm", prettyCheckbox("togglenorm_stat", HTML("<b>Toggle Normal Approximation (fixed)</b>"), TRUE)),
                     tags$div(class="norm", prettyCheckbox("togglenorm_move", "Toggle Normal Approximation (variable)", TRUE)),
                     #tags$div(class="tilt", prettyCheckbox("toggletilt_stat", HTML("<b>Toggle Tilted Distribution (fixed)</b>"), TRUE)),
                     #tags$div(class="tilt", prettyCheckbox("toggletilt_move", "Toggle Tilted Distribution (variable)", TRUE)),
                     tags$div(class="sim", prettyCheckbox("togglesim2", HTML("Toggle Simulated Distribution (visible points)"), FALSE)),
                     # textOutput("text")
                 )
               ),
               
               conditionalPanel(
                 condition = "input.tabs == 'sadd'",
                 box(width=3, tags$div(class="sim", prettyCheckbox("togglesim", HTML("<b>Toggle 'Simulated' Line</b> (all points)"), FALSE)),
                     tags$div(class="sadl", prettyCheckbox("togglesadl", HTML("<b>Toggle 'Saddlepoint' Line</b>"), TRUE)),
                     #tags$div(class="true", prettyCheckbox("toggletrue", HTML("<b>Toggle 'True' Line<b>"), TRUE))
                 )
               )
             )
    ))
)




server <- function(input, output, session) {
  
  #onclick(id="$('a.sidebar-toggle')", 
  #        shinyjs::runjs("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';"))
  
  
  n <- 0
  makeReactiveBinding('n')
  observeEvent(list(input$sidebarExpanded,input$sidebarCollapsed), { 
    n <<- n + 1
  })
  
  
  
  
  x1 <- reactiveVal()
  x2 <- reactiveVal()
  
  y <- reactiveVal()
  
  y_simulated <- reactiveVal(0)
  
  distlist1 = reactive(distfunc(input$distrString1))
  distlist2 = reactive(distfunc(input$distrString2))
  
  distr1 <- reactiveVal()
  distr2 <- reactiveVal()
  
  observeEvent(list(input$applydistr,input$distrString1, input$distrString2),{
    distr1(distlist1()$rdist)
    distr2(distlist2()$rdist)
    y(c())
    x1(c())
    x2(c())
    y_simulated(0)
    extend_points_reactive(x1, x2, y, y_simulated, 
                           input$IntensityToSimulate, 
                           distr1, distr2)
  })
  
  normalising <- reactiveVal()
  
  observeEvent(input$togglenormalising, {
    normalising(function(s){distlist1()$m(s, input$var1, input$var2) *
        distlist2()$m(s, input$var3, input$var4)})
  })
  
  observeEvent(input$IntensityToSimulate, {
    extend_points_reactive(x1, x2, y, y_simulated, 
                           input$IntensityToSimulate, 
                           distr1, distr2)
  })
  
  m = reactive({
    if (input$togglenormalising) {
      normalising()(ifelse(input$colourswitch=="true", k_prime_inverse()(input$x), input$s))
    } else {
      1
    }
  })
  
  scaled_y <- reactive({
    scale_points_y(scaled_x(), y(), 
                   ifelse(input$colourswitch=="true", k_prime_inverse()(input$x), input$s), m())
  })
  
  scaled_x = reactive({
    distlist1()$scalex(input$var1, input$var2, req(x1()))+ distlist2()$scalex(input$var3, input$var4, req(x2()))
  })
  
  
  
  
  
  
  
  
  
  get_cgfs = reactive({
    m1 =  distlist1()$m
    m2 = distlist2()$m
    m = function(sval, var1, var2, var3, var4) m1(sval, var1, var2) * m2(sval, var3, var4)
    ms = Deriv(m,"sval")
    mss = Deriv(ms, "sval")
    list("m" = m, "ms" = ms, "mss" = mss)
  })
  
  cgfs = reactive({function(s){
    cc = get_cgfs()
    var1 = input$var1; var2 = input$var2; var3 = input$var3; var4 = input$var4
    
    m = cc$m(s, var1, var2, var3, var4)
    ms = cc$ms(s, var1, var2, var3, var4)
    mss = cc$mss(s, var1, var2, var3, var4)
    
    k = log(m)
    k_prime = ms / m
    k_dblprime = (mss*m - ms^2) / m^2
    
    list("k" = k,"k_prime" = k_prime,"k_dblprime" = k_dblprime)
  }})
  
  k_prime_inverse = reactive({
    max_s = min(distlist1()$max_s(input$var1, input$var2), distlist2()$max_s(input$var3, input$var4))
    max_s = ifelse(is.finite(max_s),max_s,0)
    kpi = Vectorize(function(x){optimize(interval = c(-10,max_s), f=function(s) req(cgfs()(s)$k) - s*x)},"x")
    function(x) suppressWarnings(unlist(unname(kpi(x)[1,])))
  })
  
  s = reactive({seq(req(minmax()$min_s),req(minmax()$max_s),by=0.05)})
  xvals = reactive(seq(input$xlim[1],input$xlim[2],by=0.2))
  cgfs_input_s = reactive(cgfs()(input$s))
  cgfs_all_s = reactive(cgfs()(s()))
  sval = reactive({sval = ifelse(input$colourswitch=="true", k_prime_inverse()(input$x),input$s)
                   ifelse(is.numeric(sval), round(sval,2), sval)})
  
  
  
  output$scaled_plot_old <- renderPlot({
    req(input$s, cgfs_input_s(), input$x)
    
    par(mar=c(2, 2, 2.5, 4))
    
    plot(1,1,xlim=input$xlim, ylim=input$ylim,
         col="white",xlab="x",ylab="")
    
    title(paste("s = ", sval()), line = 1.5)
    
    #         black  yelloworange lightblue  green    yellow   darkblue  redorange purplepink
    cols = c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")[c(2,3,4,6,7)]
    pch=if(input$togglecol) 21:25 else 19
    indices = scaled_y() <= input$ylim[2]
    points(scaled_x()[indices], scaled_y()[indices],
           col=rep(cols, length.out=length(scaled_y()))[indices],
           bg=rep(cols, length.out=length(scaled_y()))[indices],
           pch=rep(pch, length.out=length(scaled_y()))[indices])
    
    lines(xvals(), y_simulated() * m() * exp((-1)*ifelse(input$colourswitch=="true", k_prime_inverse()(input$x), input$s)*xvals()),
          col = "#D55E00")
    
    if(input$togglex2){
      vline = round(ifelse(input$colourswitch=="true", input$x, cgfs()(input$s)$k_prime),2)
      abline(v = vline, col="#a5becc")
      tryCatch(mtext(paste("K'(s) = ", round(vline,2)), side=3, at = vline), error=function(e){})
    }
    
  })
  
  output$another_plot_old <- renderPlot({
    req(input$s, cgfs_input_s(), input$x)
    
    #      redorange =D55E00  0072B2      blue  lightgrey
    cols = c("black",	"#0072B2", "#E69F00","#a5becc")
    lty=if(input$togglecol) c("dashed", "solid") else rep("solid", 2)
    lwd=if(input$togglecol) c(3, 3) else c(2, 2)
    log = ifelse(input$togglelog == TRUE, "y", "")
    par(mar=c(2, 2, 2.5, 4))
    
    plot(1,1,ylab="", xlab="x",
         xlim = input$xlim, main="", type='l', log=log,
         ylim = c(max(0.00001,input$ylim2[1]), input$ylim2[2]))
    title(paste("s = ", sval()), line = 1.5)
    
    #crosshairs
    if (input$togglex){
      if(input$colourswitch=="true"){
        abline(v = input$x, col=cols[4])
        tryCatch(mtext(paste("K'(s) = ", round(input$x,2)), side=3, at = input$x), error=function(e){})
        hval = tryCatch(exp(cgfs()(k_prime_inverse()(input$x))$k - k_prime_inverse()(input$x) * input$x) /
                          suppressWarnings(sqrt(2*pi*cgfs()(k_prime_inverse()(input$x))$k_dblprime)), error=function(e){})
        abline(h=hval,col=cols[4])
        tryCatch(mtext(TeX(sprintf(r'($\frac{e^{K(s)-sx}}{\sqrt{2 \pi K''(s)}} = %f$)', round(hval,2))),
                       at=min(hval, input$ylim2[2]), side=4, padj=1), error=function(e){})
      }else{
        abline(v = cgfs_input_s()$k_prime, col=cols[4])
        tryCatch(mtext(paste("K'(s) = ", round(cgfs_input_s()$k_prime,2)), side=3, at = cgfs_input_s()$k_prime), error=function(e){})
        hval = exp(cgfs_input_s()$k - input$s * cgfs_input_s()$k_prime) / suppressWarnings(sqrt(2*pi*cgfs_input_s()$k_dblprime))
        abline(h=hval,col=cols[4])
        tryCatch(mtext(TeX(sprintf(r'($\frac{e^{K(s)-sx}}{\sqrt{2 \pi K''(s)}} = %f$)', round(hval,2))),
                       at=min(hval, input$ylim2[2]), side=4, padj=1), error=function(e){})
      }
    }
    
    #saddlepoint approximation
    x_saddle = cgfs_all_s()$k_prime
    if(input$togglesadl) lines(x_saddle, exp(cgfs_all_s()$k - s() * cgfs_all_s()$k_prime) / 
                                 suppressWarnings(sqrt(2*pi*cgfs_all_s()$k_dblprime)), 
                               col=cols[2], lwd=lwd[2], type='l', lty=lty[2])
    
    #simulated distribution
    if (input$togglesim){
      if(distlist1()$discrete & distlist2()$discrete){
        lines(table(scaled_x())/length(scaled_x()),lwd=lwd[1], col=cols[1], type='l', lty=lty[1])
      }else{
        lines(density(scaled_x(), adjust=0.8), lwd=lwd[1], col=cols[1], lty=lty[1])
      }
    }
    
    #legend
    legend(x="topright", legend = c("Simulated Distribution",  "Saddlepoint Approximation (unnormalised)"),
           col=cols, lwd=lwd, lty=lty)
    
  })
  
  output$yet_another_plot_old <- renderPlot({
    req(input$s,cgfs_input_s(), input$x)
    
    cols = c("white",	"#E69F00", "#0072B2","#a5becc", "black")
    log = ifelse(input$togglelog == TRUE, "y", "")
    lty=if(input$togglecol) c("dotted", "solid", "dashed") else rep("solid", 3)
    lwd=if(input$togglecol) c(2,4) else c(1, 2)
    par(mar=c(2, 2, 2.5, 4))
    
    plot(1,1, col=cols[1], xlim=input$xlim, type='l',
         ylim = c(max(0.00001,input$ylim2[1]), input$ylim2[2]), log=log, ylab="", xlab="x")
    title(paste("s = ", sval()), line = 1.5)
    
    
    if(input$colourswitch=="true"){
      const = 1 / suppressWarnings(sqrt(2*pi*cgfs()(k_prime_inverse()(input$x))$k_dblprime))
      kp = input$x
      if(input$togglenorm_move) lines(xvals(), dnorm(xvals(), input$x, 
                                                     suppressWarnings(sqrt(cgfs()(k_prime_inverse()(input$x))$k_dblprime))), col=cols[3], lwd=lwd[1], lty=lty[2])
      #if(input$toggletilt_move) lines(xvals(), exp(k_prime_inverse()(input$x) * xvals() - cgfs()(k_prime_inverse()(input$x))$k) *
      #                                  ddist_x(), col=cols[2], lty=lty[1], lwd=lwd[1])
    }else{
      const = 1 / suppressWarnings(sqrt(2*pi*cgfs_input_s()$k_dblprime))
      kp = cgfs_input_s()$k_prime
      if(input$togglenorm_move) lines(xvals(), dnorm(xvals(), cgfs_input_s()$k_prime, 
                                                     suppressWarnings(sqrt(cgfs_input_s()$k_dblprime))), col=cols[3], lwd=lwd[1], lty=lty[2])
      #if(input$toggletilt_move) lines(xvals(), exp(input$s * xvals() - cgfs_input_s()$k) * ddist_x(),
      #                                col=cols[2], lty=lty[1], lwd=lwd[1])
    }
    
    if(input$togglex){
      abline(v = kp, col=cols[4])
      abline(h = const, col=cols[4], lwd=1)
      mtext(paste("K'(s) = ", round(kp,2)), side=3, at = kp)
      mtext(TeX(sprintf(r'($\frac{1}{\sqrt{2\pi K''(s)}} =  %f$)', round(const,2))),
            at = min(const, input$ylim2[2]), side=4, padj=1)
    }
    
    kp2 = cgfs_all_s()$k_prime
    #if(input$toggletilt_stat) lines(kp2, exp(s() * kp2 - cgfs_all_s()$k) * ddist()(kp2), col=cols[2], lwd=lwd[2], lty=lty[1])
    if(input$togglenorm_stat) lines(kp2, 1 / suppressWarnings(sqrt(2*pi*cgfs_all_s()$k_dblprime)), col=cols[3], lwd=lwd[2], lty=lty[2])
    
    
    finite_x = scaled_x()[scaled_y() < input$ylim[2]]
    tryCatch(if (input$togglesim2){
      if(distlist1()$discrete & distlist2()$discrete){
        lines(table(finite_x)/length(finite_x), col=cols[5], type='l', lty=lty[3], lwd=lwd[1])
      }else{
        lines(density(finite_x, adjust=0.8), col=cols[5], lty=lty[3], lwd=lwd[1])
      }
    }, error=function(e){})
    
    
    
    legend(x="topright", legend = c("Normal Approximation", "Simulated Distribution"), col=cols[c(3,5)], lwd=4, lty=lty[2:3])
    
  })
  
  
  
  output$x_limit <- renderUI({sliderInput("xlim", "x Range:", min = xlimits()[1], max = xlimits()[2], value = xlimits())})
  
  xlimits = reactive({c(-25,25)})
  
  minmax = reactive({
    max  = ifelse(typeof(k_prime_inverse()(xlimits()[2]))=="list", 1, k_prime_inverse()(xlimits()[2]))
    max = min(distlist1()$max_s(input$var1, input$var2), distlist2()$max_s(input$var3, input$var4), max, xlimits()[2])
    max_x = round(cgfs()(max)$k_prime,2)
    min_x = max(distlist1()$min_x+distlist2()$min_x, xlimits()[1], round(cgfs()(-10)$k_prime,2)) 
    
    min = max(distlist1()$min_x+distlist2()$min_x, xlimits()[1])
    min_s = round(ifelse(typeof(k_prime_inverse()(min))=="list", 0, k_prime_inverse()(min)),2)
    max  = ifelse(typeof(k_prime_inverse()(xlimits()[2]))=="list", 1, k_prime_inverse()(xlimits()[2]))
    max_s = round(min(distlist1()$max_s(input$var1, input$var2),distlist2()$max_s(input$var3, input$var4), max),2)
    
    list("max_s" = max_s, "max_x" = max_x, "min_x" = min_x, "min_s" = min_s)
  })
  
  output$xLoopSlider <- renderUI({
    conditionalPanel(condition = "input.colourswitch == 'true'",
                     sliderInput("x", "Select K'(s):",
                                 min = req(minmax()$min_x), max=req(minmax()$max_x),
                                 value = 1, step = 0.05, width="100%",
                                 animate = animationOptions(interval = 300, loop = TRUE)))
  })
  
  
  output$sLoopSlider <- renderUI({
    conditionalPanel(condition = "input.colourswitch == 'false'",
                     sliderInput("s", "Select s:",
                                 min = req(minmax()$min_s), max=req(minmax()$max_s),
                                 value = 0, step = 0.05,width="100%",
                                 animate = animationOptions(interval = 300, loop = TRUE)))
  })
  
  
  
  observeEvent(input$x, if(input$colourswitch=="true"){
    updateNumericInput(session, inputId = "s", value = k_prime_inverse()(input$x))})
  
  observeEvent(list(input$var1, input$var2),
               if(input$colourswitch=="false"){updateNumericInput(session, inputId = "s", value = input$s)})
  
  observeEvent(input$s, if(input$colourswitch=="false"){
    updateNumericInput(session, inputId = "x", value = cgfs()(input$s)$k_prime)})
  
  observeEvent(list(input$var1, input$var2),
               if(input$colourswitch=="true"){updateNumericInput(session, inputId = "x", value = input$x)}) #look. it works. and breaks without it. so.
  
  
  
  output$name_first_slider <- renderUI(HTML(paste("&emsp;<b>",distlist1()$slidernames[1])))
  
  output$name_second_slider <- renderUI(HTML(paste("&emsp;<b>",distlist1()$slidernames[2])))
  
  output$name_third_slider <- renderUI(HTML(paste("&emsp;<b>",distlist2()$slidernames[1])))
  
  output$name_fourth_slider <- renderUI(HTML(paste("&emsp;<b>",distlist2()$slidernames[2])))
  
  output$first_slider <- renderUI({
    sliderInput("var1", label= NULL, step = distlist1()$step[1],
                min = distlist1()$firstminmax[1], max = distlist1()$firstminmax[2],
                value = distlist1()$start[1], animate =
                  animationOptions(interval = 1000, loop = TRUE))
  })
  
  output$second_slider <- renderUI({
    conditionalPanel(condition = distlist1()$condition,
                     sliderInput("var2", label= NULL, step=distlist1()$step[2],
                                 min = distlist1()$secondminmax(input$var1, input$var2)[1],
                                 max = distlist1()$secondminmax(input$var1, input$var2)[2],
                                 value = distlist1()$start[2],animate =
                                   animationOptions(interval = 1000, loop = TRUE))
    )
  })
  
  output$third_slider <- renderUI({
    sliderInput("var3", label= NULL, step = distlist2()$step[1],
                min = distlist2()$firstminmax[1], max = distlist2()$firstminmax[2],
                value = distlist2()$start[1], animate =
                  animationOptions(interval = 1000, loop = TRUE))
  })
  
  output$fourth_slider <- renderUI({
    conditionalPanel(condition = distlist2()$condition,
                     sliderInput("var4", label= NULL, step=distlist2()$step[2],
                                 min = distlist2()$secondminmax(input$var3, input$var4)[1],
                                 max = distlist2()$secondminmax(input$var3, input$var4)[2],
                                 value = distlist2()$start[2],animate =
                                   animationOptions(interval = 1000, loop = TRUE))
    )
  })

  output$changesidebartoggleicon <- renderUI({
    if(n > 1){
      tags$head(tags$style(HTML('
      .main-header .sidebar-toggle:before {
        content: \' \\f0c9 \';}')))
    }
  })
  
  output$slidercols <- renderUI({
    colour <- "#E69F00"
    tagList(
      tags$style(sprintf(".irs-bar-edge, .irs-bar, .irs-single, .irs-from, .irs-to {background: %s !important;}", colour)),
      tags$style(sprintf(".irs-bar {border-top: %s !important;}", colour)),
      tags$style(sprintf(".irs-bar {border-bottom: %s !important;}", colour))
    )
  })
  
}

shinyApp(ui, server)