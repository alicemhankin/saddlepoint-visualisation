library(shinydashboard)
library(shiny)
library(latex2exp)
library(shinyWidgets)
library(glue)
library(Deriv)
options(shiny.useragg = TRUE) 



#########################################################################################################################################################
# This application was written by Alice Hankin, adapted from code by Jesse Goodman and Rachel Fewster                                                   #
# More information about this app, along with my thesis about this project can be found at https://github.com/alicemhankin/saddlepoint-visualisation    #
# A start guide can be found at https://raw.githack.com/alicemhankin/saddlepoint-visualisation/main/user-guide.html                                     #
#########################################################################################################################################################




#function returns info pertaining to the selected distribution
distfunc = function(distrString, min_x, max_s){

  names = c("ddist", "rdist", "m", "firstminmax", "secondminmax", "step", "max_s", "slidernames", 
            "scalex", "discrete", "condition", "min_x")
  mylist <- sapply(names,function(x) NULL)
  
  #The members of mylist are (in order):
  #true distribution,  initial sampling distribution,  moment generating function,
  #minimum/maximum of first parameter, minimum/maximum of second parameter,  step size for parameters, 
  #maximum value s can take, names of parameter sliders, values from rdist adapted for chosen parameters,
  #is sampling distribution discrete?, are there two parameters?,  minimum value of x
  
  #note that the initially sampled points from rdist are fixed (until we resample or change distribution)
  #but scalex describes how to adapt these values so that they follow the correct distribution with the correct 
  #parameters. For example, for Poisson, rdist is runif(n) which are not Poisson distributed. scalex translates
  #the random uniform variables x into Poisson(a) variables via qpois(x, a).
  
  #p.s. many of these parameters are included just to avoid minor problems (especially purely aesthetic
  #ones). I didn't bother trying to decrease the number of inputs here or make it super easy to add new distributions!
  # -- from the 'bivariate' app we see that all that is /really/ needed inputted here is m and rdist!!
  
  #default step size
  def_step = 0.05
  
  mylist$min_x = min_x
  
  switch(distrString,
         
         "Normal" = {mylist$start = c(0,1);
         mylist$ddist = function(x, var1, var2) dnorm(x, var1, var2);
         mylist$rdist = function(n) rnorm(n, mean=mylist$start[1],
                                          sd = mylist$start[2]);
         mylist$m = function(s, var1, var2) exp(s*var1+0.5*s^2*var2^2);
         mylist$firstminmax = c(-10,10);
         mylist$secondminmax = function(var1, var2) c(def_step,3.4);
         mylist$step = rep(def_step, 2);
         mylist$max_s = function(var1, var2) max_s;
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
         mylist$ddist = function(x, var1, var2) ifelse(x > var2+1 | x < -1 ,NA,suppressWarnings(factorial(var2) * var1^x * (1-var1)^(var2-x)/(factorial(x)*factorial(var2-x))));
         mylist$rdist = function(n) runif(n);
         mylist$m = function(s, var1, var2) (1-var1+var1*exp(s))^var2;
         mylist$firstminmax = c(def_step,1-def_step);
         mylist$secondminmax = function(var1, var2) c(1,10);
         mylist$step = c(def_step,1);
         mylist$max_s = function(var1, var2) log((-var2* var1 + 0.05 *var1 + var2 - 0.05)/(0.05*var1));
         mylist$slidernames = c("Probability of a success:", "Number of Trials:");
         mylist$scalex = function(new1,new2,x) suppressWarnings(qbinom(x, new2, new1));
         mylist$discrete = T
         mylist$min_x = 0.05 # this 0.05 is the same as in max_s
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
         mylist$max_s = function(var1, var2) max_s;
         mylist$slidernames = c("Minimum:", "Maximum:");
         mylist$scalex = function(new1,new2,x) suppressWarnings(qunif(punif(x,mylist$start[1],mylist$start[2]),new1,new2));
         mylist$discrete = F
         }
         
  )
  
  mylist$condition = ifelse(mylist$slidernames[2]=="", "false", "true")
  
  mylist
}

#scaling y coordinates
scale_points_y = function(x_points, y_points, s, m) y_points*m*exp(-s*x_points)

#defines how to create new points
extend_points_reactive = function(x, y, y_simulated, simulate_to_y, distr) {
  # x is a vector of numbers drawn from distr
  # y is a vector of numbers drawn from a rate 1 Poisson point process on (0,infty)
  # y_simulated is the value for which the P.p.p. has been simulated 
  if (simulate_to_y > y_simulated()) {
    new_count = rpois(1, simulate_to_y - y_simulated())
    new_x = distr()(new_count)
    new_y = runif(new_count, min = y_simulated(), max = simulate_to_y)
    x(c(x(), new_x))
    y(c(y(), new_y))
    y_simulated(simulate_to_y)
  }
}


ui <- dashboardPage(
  
  #navigation bar
  dashboardHeader(
    title = "A Visual Representation of the Saddlepoint Approximation",
    titleWidth = 600),
  
  #collapsible sidebar
  dashboardSidebar(
    
    width=300,
    
    #dropdown menu for choosing the distribution of the x coordinates
    selectInput("distrString", label = "Distribution for x:", 
                choices = c("Normal", "Exponential","Gamma", "Chi-Squared", "Poisson",
                            "Geometric", "Negative Binomial", "Binomial", "Uniform"), width="100%"),
    
    htmlOutput("name_first_slider"),
    
    uiOutput("first_slider"),
    
    htmlOutput("name_second_slider"),
    
    uiOutput("second_slider"),
    
    #colour of the re-simulate button matches the navigation bar
    actionButton("applydistr", "Re-Simulate", icon("refresh", lib="glyphicon"),
                 style="color: #FFFFFF; background-color: var(--nav-color)"),
    
    #horizontal line
    hr(),
    
    #buttons which allow selection of mode -- does the user want to choose s or K'(s)?
    radioGroupButtons(
      inputId = "colourswitch",
      choices = c( "Choose s"="false", "Choose K'(s)"="true"),
      justified=TRUE, label=NULL
    ),
    
    uiOutput("xLoopSlider"),
    
    uiOutput("sLoopSlider"),
    
    uiOutput("x_limit"),
    
    #toggle for colourblind mode
    prettyCheckbox("togglecol", "Toggle Colourblind Mode", FALSE)
    
  ),
  
  #main panel
  dashboardBody(
    
    #removing scrollbar
    tags$head(tags$style(".body {overflow-y: hidden;}")),
    
    #changing spacing between widgets
    tags$style(".form-group {margin-bottom: 9px;}"),
    
    #changing colour of radio button labels (toggle graph lines on/off)
    tags$style(".tilt {color: #E69F00;}"),
    tags$style(".norm {color: #0072B2;}"),
    tags$style(".sim {color: #000000;}"),
    tags$style(".sadl {color: #0072B2;}"),
    tags$style(".true {color: #E69F00;}"),
    
    #changing colour of navigation bar
    uiOutput("navcolor"),
    tags$head(tags$style(HTML('
      .skin-blue .main-header .logo {background-color: var(--nav-color);}
      .skin-blue .main-header .logo:hover {background-color: var(--nav-color);}
      .skin-blue .main-header .navbar .sidebar-toggle:hover{background-color: var(--nav-color);}
      .skin-blue .main-header .navbar {background-color: var(--nav-color);}'))),
    
    #changing colour of active tab indicator to match navigation bar
    tags$style(".nav-tabs-custom .nav-tabs li.active {border-top-color: var(--nav-color);}"),
    
    #setting background colour
    setBackgroundColor("#efefef", shinydashboard = T),
    
    uiOutput("slidercols"),
    
    uiOutput("changesidebartoggleicon"),
    
    #upper row in main panel
    fluidRow(
      
      #scatterplot box
      box(width=9, plotOutput("scatterplot", height = "38vh"), height = "41vh"),
      
      #controls for scatterplot
      box(width=3,
          sliderInput("ylim", "y Range:",min = 0, max = 1000,value = c(0,100)),
          sliderInput("IntensityToSimulate", "Intensity to simulate:",
                      min = 0, max = 100000, value = 100),
          fluidRow(column(6,prettyCheckbox("togglex2", "Toggle Crosshairs", TRUE)),
                   column(6,prettyCheckbox("togglenormalising", "Normalise", TRUE))))),
    
    #lower row in main panel
    fluidRow(
      
      #tab box for saddlepoint approximation plot and tilted density plot
      tabBox(width=9,id="tabs",
             tabPanel("Saddlepoint Approximation", value="sadd", plotOutput("saddlepoint", height = "39vh")),
             tabPanel("Tilted Distribution",  value="tilt",plotOutput("tilted",height = "39vh")), height = "42vh"),
      
      #joint controls for saddlepoint and tilted plots
      box(width=3, 
          sliderInput("ylim2", "y Range:", min = 0, max = 2, value = c(0,0.5), step=0.01),
          fluidRow(column(6,prettyCheckbox("togglex", "Toggle Crosshairs", TRUE)),
                   column(6,prettyCheckbox("togglelog", "Toggle Log", FALSE)))
      ),
      
      #controls for tilted density plot only appear if that tab is selected
      conditionalPanel(
        condition = "input.tabs == 'tilt'",
        box(width=3, 
            tags$div(class="tilt", prettyCheckbox("toggletilt_move", "Tilted Distribution (selected s)", TRUE)),
            tags$div(class="tilt", prettyCheckbox("toggletilt_stat", HTML("<b>Mean Tilted Distribution (all s)</b>"), TRUE)),
            tags$div(class="norm", prettyCheckbox("togglenorm_move", "Normal Approximation (selected s)", TRUE)),
            tags$div(class="norm", prettyCheckbox("togglenorm_stat", HTML("<b>Mean Normal Approximation (all s)</b>"), TRUE)),
            tags$div(class="sim", prettyCheckbox("togglesim2", HTML("Simulated Distribution (visible points)"), FALSE))
        )
      ),
      
      #controls for saddlepoint approximation plot only appear if that tab is selected
      conditionalPanel(
        condition = "input.tabs == 'sadd'",
        box(width=3, 
            tags$div(class="true", prettyCheckbox("toggletrue", HTML("<b>Toggle 'True' Line<b>"), TRUE)),
            tags$div(class="sadl", prettyCheckbox("togglesadl", HTML("<b>Toggle 'Saddlepoint' Line</b>"), TRUE)),
            tags$div(class="sim", prettyCheckbox("togglesim", HTML("<b>Toggle 'Simulated' Line</b> (all points)"), FALSE))
        )
      )
    )
  )
)


server <- function(input, output, session) {
  
  #initializes vectors for x and y coordinates of points, and sampling distribution
  x <- reactiveVal()
  y <- reactiveVal()
  y_simulated <- reactiveVal(0)
  distr <- reactiveVal()
  
  #gets all information about selected distribution
  distlist = reactive(distfunc(input$distrString, req(minmax_x()[1]), req(minmax_s()[2])))

  #updates x values, y values, and sampling distribution
  observeEvent(list(input$applydistr,input$distrString),{
    distr(distlist()$rdist)
    y(c())
    x(c())
    y_simulated(0)
    extend_points_reactive(x, y, y_simulated, 
                           input$IntensityToSimulate, 
                           distr)
  })
  
  #updates mgf when the normalising toggle is toggled
  normalising <- reactiveVal()
  observeEvent(input$togglenormalising, {
    normalising(function(s){distlist()$m(s, input$var1, input$var2)})
  })
  
  #adds more points if intensity to simulate is increased
  observeEvent(input$IntensityToSimulate, {
    extend_points_reactive(x, y, y_simulated, 
                           input$IntensityToSimulate, 
                           distr)
  })
  
  #returns factor to multiply y by -- this is the mgf if normalising is on, and 1 if it is off
  m = reactive({
    if (input$togglenormalising) {
      normalising()(ifelse(input$colourswitch=="true", k_prime_inverse()(input$x), input$s))
    } else {
      1
    }
  })
  
  #returns final y coordinates (normalised iff toggle is turned on)
  scaled_y <- reactive({
    scale_points_y(scaled_x(), y(), ifelse(input$colourswitch=="true", k_prime_inverse()(input$x), input$s), m())
  })
  
  #returns final x coordinates
  scaled_x = reactive({
    distlist()$scalex(input$var1, input$var2, req(x()))
  })
  
  #finds derivatives of moment generating function
  get_cgfs = reactive({
    m= distlist()$m
    ms = Deriv(m,"s")
    mss = Deriv(ms, "s")
    list("m" = m, "ms" = ms, "mss" = mss)
  })
  
  #finds derivatives of cumulant generating function
  cgfs = reactive({function(s){
    c = get_cgfs()
    var1 = input$var1; var2 = input$var2
    
    m = c$m(s, var1, var2)
    ms = c$ms(s, var1, var2)
    mss = c$mss(s, var1, var2)
    
    k = log(m)
    k_prime = ms / m
    k_dblprime = (mss*m - ms^2) / m^2
    
    #agghhhh the mgf of uniform dist is piecewise. I know this is the worst but ¯\_(ツ)_/¯
    if(input$distrString=="Uniform"){k[s==0] = 0; k_prime[s==0] = 0.5*(var2+var1); k_dblprime[s==0] = (var1-var2)^2/12} 
    
    list("k" = k,"k_prime" = k_prime,"k_dblprime" = k_dblprime)
  }})
  
  #a function which is the inverse of K'(s)
  k_prime_inverse = reactive({
    max_s = distlist()$max_s(input$var1, input$var2)
    max_s = ifelse(is.finite(max_s),max_s,0)
    kpi = Vectorize(function(x){optimize(interval = c(minmax_s()[1],max_s), f=function(s) req(cgfs()(s)$k) - s*x)},"x")
    function(x) suppressWarnings(unlist(unname(kpi(x)[1,])))
  })
  
  #defines values that are needed several times in plots, so they don't have to be run multiple times
  s = reactive({seq(req(minmax()$min_s),req(minmax()$max_s),by=0.01)})
  xvals = reactive(seq(input$xlim[1],input$xlim[2],by=0.01))
  ddist = reactive(function(x) distlist()$ddist(x, input$var1, input$var2))
  ddist_x = reactive(ddist()(xvals()))
  cgfs_input_s = reactive(cgfs()(input$s))
  cgfs_all_s = reactive(cgfs()(s()))
  
  
  
  #draws scatterplot
  output$scatterplot <- renderPlot({
    req(input$s, cgfs_input_s(), input$x)
    
    #current s value
    sval = ifelse(input$colourswitch=="true", k_prime_inverse()(input$x),input$s)
    sval = ifelse(is.numeric(sval), round(sval,2), sval)
    
    #draws axes and title
    par(mar=c(2, 2, 2.5, 4))
    plot(1,1,xlim=input$xlim, ylim=input$ylim, main=paste("s = ", sval), col="white",xlab="x",ylab="")
    
    #selects colours and shapes
    #         black  yellow-orange lightblue  green    yellow   darkblue  red-orange purple-pink
    cols = c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")[c(2,3,4,6,7)]
    pch=if(input$togglecol) 21:25 else 19
    
    #which points should be drawn?
    indices = tryCatch((scaled_y() <= input$ylim[2]), error=function(e){})
    
    #draws points
    tryCatch(points(scaled_x()[indices], scaled_y()[indices],
                    col=rep(cols, length.out=length(scaled_y()))[indices],
                    bg=rep(cols, length.out=length(scaled_y()))[indices],
                    pch=rep(pch, length.out=length(scaled_y()))[indices]), error=function(e){})
    
    #maximum y value that points can have (red line)
    tryCatch(lines(xvals(), y_simulated() * m() * exp((-1)*ifelse(input$colourswitch=="true", 
                                                                  k_prime_inverse()(input$x), input$s)*xvals()),
                   col = "#D55E00"), error=function(e){})
    
    #draws crosshairs
    if(input$togglex2){
      vline = round(ifelse(input$colourswitch=="true", input$x, cgfs_input_s()$k_prime),2)
      abline(v = vline, col="#a5becc")
      tryCatch(mtext(paste("K'(s) = ", round(vline,2)), side=3, at = vline), error=function(e){})
    }
    
  })
  
  #draws saddlepoint approximation plot
  output$saddlepoint <- renderPlot({
    req(input$s, cgfs_input_s(), input$x)
    
    #selects colours, line weights, line types, and y scale
    #                 dark-blue  yellow-orange  light-grey-blue
    cols = c("black",	"#0072B2", "#E69F00",     "#a5becc")
    lty=if(input$togglecol) c("dashed", "solid", "dotted") else rep("solid", 3)
    lwd=if(input$togglecol) c(3, 3, 3) else c(2, 2, 2)
    log = ifelse(input$togglelog == TRUE, "y", "")
    
    #draws plot axes
    par(mar=c(2, 2, 1, 4))
    plot(1,1,ylab="", xlab="x",
         xlim = input$xlim, main="", type='l', log=log,
         ylim = c(max(0.00001,input$ylim2[1]), input$ylim2[2]))
    
    #draws crosshairs
    if (input$togglex){
      if(input$colourswitch=="true"){
        abline(v = input$x, col=cols[4])
        tryCatch(mtext(paste("K'(s) = ", round(input$x,2)), side=3, at = input$x), error=function(e){})
        hval = tryCatch(exp(cgfs()(k_prime_inverse()(input$x))$k - k_prime_inverse()(input$x) * input$x) /
                          suppressWarnings(sqrt(2*pi*cgfs()(k_prime_inverse()(input$x))$k_dblprime)), error=function(e){})
      }else{
        abline(v = cgfs_input_s()$k_prime, col=cols[4])
        tryCatch(mtext(paste("K'(s) = ", round(cgfs_input_s()$k_prime,2)), side=3, at = cgfs_input_s()$k_prime), error=function(e){})
        hval = exp(cgfs_input_s()$k - input$s * cgfs_input_s()$k_prime) / suppressWarnings(sqrt(2*pi*cgfs_input_s()$k_dblprime))
      }
      abline(h=hval,col=cols[4])
      tryCatch(mtext(TeX(sprintf(r'($\frac{e^{K(s)-sx}}{\sqrt{2 \pi K''(s)}} = %f$)', round(hval,2))),
                     at=min(hval, 0.8*input$ylim2[2]), side=4, padj=1), error=function(e){})
    }
    
    #true distribution
    if(input$toggletrue) lines(xvals(),  ddist_x(),
                               col=cols[3], lty=lty[3], lwd=lwd[3])
    
    #saddlepoint approximation
    x_saddle = cgfs_all_s()$k_prime
    if(input$togglesadl) lines(x_saddle, exp(cgfs_all_s()$k - s() * cgfs_all_s()$k_prime) / 
                                 suppressWarnings(sqrt(2*pi*cgfs_all_s()$k_dblprime)), 
                               col=cols[2], lwd=lwd[2], type='l', lty=lty[2])
    
    #simulated distribution
    tryCatch(if (input$togglesim){
      if(distlist()$discrete){
        lines(table(scaled_x())/length(scaled_x()),lwd=lwd[1], col=cols[1], type='l', lty=lty[1])
      }else{
        lines(density(scaled_x(), adjust=0.8), lwd=lwd[1], col=cols[1], lty=lty[1])
      }
    }, error=function(e){})
    
    #legend
    legend(x="topright", legend = c("Simulated Distribution",  "Saddlepoint Approximation (unnormalised)",
                                    "True Distribution"), col=cols, lwd=lwd, lty=lty)
    
  })
  
  #draws tilted density plot
  output$tilted <- renderPlot({
    req(input$s,cgfs_input_s(), input$x)
    
    #selects colours, line weights, line types, and y scale
    #                 yellow-orange  dark-blue  light-grey-blue
    cols = c("white",	"#E69F00",     "#0072B2", "#a5becc",       "black")
    log = ifelse(input$togglelog == TRUE, "y", "")
    lty=if(input$togglecol) c("dotted", "solid", "dashed") else rep("solid", 3)
    lwd=if(input$togglecol) c(2,4) else c(1, 2)
    
    #initializes plot
    par(mar=c(2, 2, 1, 4))
    plot(1,1, col=cols[1], xlim=input$xlim, type='l',
         ylim = c(max(0.00001,input$ylim2[1]), input$ylim2[2]), log=log, ylab="", xlab="x")
    
    #tilted density plot and normal approximation to the tilted density
    if(input$colourswitch=="true"){
      const = 1 / suppressWarnings(sqrt(2*pi*cgfs()(k_prime_inverse()(input$x))$k_dblprime))
      kp = input$x
      if(input$togglenorm_move) lines(xvals(), dnorm(xvals(), input$x, 
                                                     suppressWarnings(sqrt(cgfs()(k_prime_inverse()(input$x))$k_dblprime))), col=cols[3], lwd=lwd[1], lty=lty[2])
      if(input$toggletilt_move) lines(xvals(), exp(k_prime_inverse()(input$x) * xvals() - cgfs()(k_prime_inverse()(input$x))$k) *
                                        ddist_x(), col=cols[2], lty=lty[1], lwd=lwd[1])
    }else{
      const = 1 / suppressWarnings(sqrt(2*pi*cgfs_input_s()$k_dblprime))
      kp = cgfs_input_s()$k_prime
      if(input$togglenorm_move) lines(xvals(), dnorm(xvals(), cgfs_input_s()$k_prime, 
                                                     suppressWarnings(sqrt(cgfs_input_s()$k_dblprime))), col=cols[3], lwd=lwd[1], lty=lty[2])
      if(input$toggletilt_move) lines(xvals(), exp(input$s * xvals() - cgfs_input_s()$k) * ddist_x(),
                                      col=cols[2], lty=lty[1], lwd=lwd[1])
    }
    
    #crosshairs
    if(input$togglex){
      abline(v = kp, col=cols[4])
      abline(h = const, col=cols[4], lwd=1)
      mtext(paste("K'(s) = ", round(kp,2)), side=3, at = kp)
      mtext(TeX(sprintf(r'($\frac{1}{\sqrt{2\pi K''(s)}} =  %f$)', round(const,2))),
            at = min(const, 0.8*input$ylim2[2]), side=4, padj=1)
    }
    
    # mean tilted and mean normal distributions for all s
    kp2 = cgfs_all_s()$k_prime
    if(input$toggletilt_stat) lines(kp2, exp(s() * kp2 - cgfs_all_s()$k) * ddist()(kp2), col=cols[2],
                                    lwd=lwd[2], lty=lty[1])
    if(input$togglenorm_stat) lines(kp2, 1 / suppressWarnings(sqrt(2*pi*cgfs_all_s()$k_dblprime)),
                                    col=cols[3], lwd=lwd[2], lty=lty[2])
    
    #simulated density
    finite_x = scaled_x()[scaled_y() < input$ylim[2]]
    tryCatch(if (input$togglesim2){
      if(distlist()$discrete){
        lines(table(finite_x)/length(finite_x), col=cols[5], type='l', lty=lty[3], lwd=lwd[1])
      }else{
        lines(density(finite_x, adjust=0.8), col=cols[5], lty=lty[3], lwd=lwd[1])
      }
    }, error=function(e){})
    
    #legend
    legend(x="topright", legend = c("Tilted distribution",
                                    "Normal Approximation", "Simulated Distribution"), col=cols[c(2,3,5)], lwd=4, lty=lty)
    
  })
  
  
  
  #minimum and maximum possible allowed values of x=K'(s) and s
  minmax_x = reactive(c(-10,15))
  minmax_s = reactive(c(-10,10))
  
  #x limit slider
  output$x_limit <- renderUI({sliderInput("xlim", "x Range:", min = minmax_x()[1], max = minmax_x()[2], value = c(-10,10))})
  
  #minimum and maximum allowed values of s and K'(s) GIVEN the distribution and parameters
  minmax = reactive({
    max  = ifelse(typeof(k_prime_inverse()(minmax_x()[2]))=="list", 1, k_prime_inverse()(minmax_x()[2]))
    max = min(distlist()$max_s(input$var1, input$var2), max, minmax_s()[2])
    max_x = round(cgfs()(max)$k_prime,2)
    min_x = max(distlist()$min_x, minmax_x()[1], round(cgfs()(minmax_s()[1])$k_prime,2)) 
    
    min = max(distlist()$min_x, minmax_x()[1])
    min_s = round(ifelse(typeof(k_prime_inverse()(min))=="list", 0, k_prime_inverse()(min)),2)
    max  = ifelse(typeof(k_prime_inverse()(minmax_x()[2]))=="list", 1, k_prime_inverse()(minmax_x()[2]))
    max_s = round(min(distlist()$max_s(input$var1, input$var2), max),2)
    
    list("max_s" = max_s, "max_x" = max_x, "min_x" = min_x, "min_s" = min_s)
  })
  
  #input for K'(s)
  output$xLoopSlider <- renderUI({
    conditionalPanel(condition = "input.colourswitch == 'true'",
                     sliderInput("x", "Select K'(s):",
                                 min = req(minmax()$min_x), max=req(minmax()$max_x),
                                 value = 1, step = 0.05,
                                 animate = animationOptions(interval = 300, loop = TRUE)))
  })
  
  #input for s
  output$sLoopSlider <- renderUI({
    conditionalPanel(condition = "input.colourswitch == 'false'",
                     sliderInput("s", "Select s:",
                                 min = req(minmax()$min_s), max=req(minmax()$max_s),
                                 value = 0, step = 0.05,
                                 animate = animationOptions(interval = 300, loop = TRUE)))
  })
  

  
  #these ensure that the K'(s) and s sliders correspond to each other
  observeEvent(input$x, if(input$colourswitch=="true"){
    updateNumericInput(session, inputId = "s", value = k_prime_inverse()(input$x))})
  
  observeEvent(list(input$var1, input$var2),
               if(input$colourswitch=="false"){updateNumericInput(session, inputId = "s", value = input$s)})
  
  observeEvent(input$s, if(input$colourswitch=="false"){
    updateNumericInput(session, inputId = "x", value = cgfs()(input$s)$k_prime)})
  
  observeEvent(list(input$var1, input$var2),
               if(input$colourswitch=="true"){updateNumericInput(session, inputId = "x", value = input$x)}) #look. it works. and breaks without it. 
  
  
  #label for first parameter input
  output$name_first_slider <- renderUI(HTML(paste("&emsp;<b>",distlist()$slidernames[1])))
  
  #label for second parameter input
  output$name_second_slider <- renderUI(HTML(paste("&emsp;<b>",distlist()$slidernames[2])))
  
  #input for first parameter
  output$first_slider <- renderUI({
    sliderInput("var1", label= NULL, step = distlist()$step[1],
                min = distlist()$firstminmax[1], max = distlist()$firstminmax[2],
                value = distlist()$start[1], animate =
                  animationOptions(interval = 500, loop = TRUE))
  })
  
  #input for second parameter
  output$second_slider <- renderUI({
    conditionalPanel(condition = distlist()$condition,
                     sliderInput("var2", label= NULL, step=distlist()$step[2],
                                 min = distlist()$secondminmax(input$var1, input$var2)[1],
                                 max = distlist()$secondminmax(input$var1, input$var2)[2],
                                 value = distlist()$start[2],animate =
                                   animationOptions(interval = 500, loop = TRUE))
    )
  })
  
  #changes colour of navbar
  output$navcolor <- renderUI({
    tags$style(HTML(glue(
      ":root {--nav-color: @{'#D55E00'}@}",
      .open = "@{", .close = "}@"
    )))
  })
  
  #changes colour of sliders
  output$slidercols <- renderUI({
    colour <- "#D55E00"
    tagList(
      tags$style(sprintf(".irs-bar-edge, .irs-bar, .irs-single, .irs-from, .irs-to {background: %s !important;}", colour)),
      tags$style(sprintf(".irs-bar {border-top: %s !important;}", colour)),
      tags$style(sprintf(".irs-bar {border-bottom: %s !important;}", colour))
    )
  })
  
  
  #defines counter for how many times the sidebar has been opened
  n <- 0
  makeReactiveBinding('n')
  observeEvent(list(input$sidebarExpanded,input$sidebarCollapsed), { 
    n <<- n + 1
  })
  
  #changes sidebar toggle icon based on how many times the sidebar has been opened
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
  
  
  #defines modal dialog
  query_modal <- modalDialog(
    title = "Have you read the start guide?", 
    footer = modalButton("I already know what I'm doing"),
    easyClose=T,
    tagList(a("Click here to learn about to use this application", target="_blank", 
              href="https://rawcdn.githack.com/alicemhankin/saddlepoint-visualisation/d1e0456747bf680a1506b686c4602583834ab86b/user-guide.html"))
    
  )
  
  #show the modal dialog on start up
  showModal(query_modal)
  
}

shinyApp(ui, server)