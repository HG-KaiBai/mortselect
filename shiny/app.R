# Native Shiny companion. From the repository root: shiny::runApp("shiny")
library(shiny)
library(ggplot2)
library(jsonlite)

data_dir <- normalizePath("../docs/data", mustWork = FALSE)
if (!dir.exists(data_dir)) data_dir <- normalizePath("../mortselect/docs/data", mustWork = TRUE)
manifest <- read_json(file.path(data_dir, "manifest.json"), simplifyVector = FALSE)
methods <- setNames(vapply(manifest$methods, `[[`, "", "id"), vapply(manifest$methods, `[[`, "", "label"))
forecasts <- setNames(lapply(unname(methods), function(id) read_json(file.path(data_dir, paste0(id, ".json")), simplifyVector = TRUE)), unname(methods))
pretty_population <- function(x) {
  x <- sub("^-- *", "", x)
  x <- sub("England & Wales Total Population", "England & Wales", x, fixed = TRUE)
  x <- sub("France total population", "France", x, fixed = TRUE)
  sub("NZ total population", "New Zealand", x, fixed = TRUE)
}
addResourcePath("explorer-assets", dirname(data_dir))
theme_css <- "
body{background:#f7f8f4;color:#183c38;font-family:'Segoe UI',Arial,sans-serif} .container-fluid{max-width:1400px;padding:0 36px}
.app-header{padding:24px 0;border-bottom:1px solid #dce3dc;display:flex;justify-content:space-between;align-items:center}.app-brand{font-size:26px;font-weight:650;letter-spacing:-1px}
a{color:#247466}.hero-title{font:48px/1.1 Georgia,serif;margin:36px 0 18px}.hero-title em{color:#247466}.intro{color:#697774;max-width:690px;line-height:1.8;margin-bottom:30px}
.well{background:#eef2ec;border:1px solid #dbe2d9;box-shadow:none;border-radius:12px;padding:23px}.well h3{font-size:17px;margin:0 0 22px}.control-label{font-size:12px}.form-control,.selectize-input{border-color:#cdd8ce;box-shadow:none;border-radius:6px}.help-block{font-size:11px;color:#697774}.radio,.checkbox{font-size:12px}.irs--shiny .irs-bar,.irs--shiny .irs-single{background:#247466;border-color:#247466}
.nav-pills{background:#e9eee6;border-radius:8px;padding:5px;margin-bottom:18px;display:inline-flex}.nav-pills>li>a{padding:10px 22px;border-radius:5px;color:#697774}.nav-pills>li.active>a,.nav-pills>li.active>a:hover{background:white;color:#183c38;box-shadow:0 1px 4px #183c3812}
.plot-panel{background:white;border:1px solid #dce3dc;border-radius:12px;padding:24px}.plot-panel h3{margin:0 0 10px;font-size:22px}.plot-subtitle{color:#697774;font-size:12px}.actions{display:flex;gap:12px;align-items:center;margin:18px 0}.btn-default{background:#183c38;color:white;border:0;border-radius:6px;padding:9px 14px}.note{border-top:1px solid #dce3dc;margin-top:30px;padding:24px 0;color:#697774;font-size:12px;line-height:1.8}.shiny-output-error-validation{color:#697774;font-size:15px}.table{font-size:12px}.tabbable{overflow:auto}
@media(max-width:700px){.container-fluid{padding:0 16px}.hero-title{font-size:36px}.plot-panel{padding:14px}.nav-pills>li>a{padding:9px 12px}.app-header{gap:12px;font-size:12px}}
"
ui <- fluidPage(
  tags$head(tags$title("mortselect — Mortality Forecast Explorer"), tags$link(rel="icon",href="explorer-assets/favicon.svg"), tags$style(HTML(theme_css))),
  div(class="app-header",div(class="app-brand","mortselect."),div(a("Public explorer ↗",href="https://hg-kaibai.github.io/mortselect/",target="_blank")," · ",a("R package ↗",href="https://github.com/HG-KaiBai/mortselect",target="_blank"))),
  h1(class="hero-title","A clearer view of ",em("mortality forecasts.")),
  p(class="intro","Explore saved research forecasts across ages, populations, and five selection and ensemble methods. Forecast horizon: 2011–2025."),
  sidebarLayout(
    sidebarPanel(width=3,h3("Build your view"),selectInput("method","01  Forecast method",choices=methods,selected="DSA"),
      uiOutput("population_control"),uiOutput("sex_control"),uiOutput("dimension_control"),checkboxInput("log","Log scale",FALSE),helpText(textOutput("coverage"))),
    mainPanel(width=9,
      tabsetPanel(id="mode",type="pills",tabPanel("By age",value="age"),tabPanel("By population",value="population"),tabPanel("By year",value="year")),
      div(class="plot-panel",h3(textOutput("title")),p(class="plot-subtitle",textOutput("subtitle")),plotOutput("forecast_plot",height="430px")),
      div(class="actions",downloadButton("download_csv","Download CSV"),downloadButton("download_plot","Download plot"),checkboxInput("table","View data",FALSE)),
      conditionalPanel("input.table",tableOutput("values"))
    )
  ),
  div(class="note",strong("A note on the data. "),"These are saved point forecasts, not observed mortality rates or newly fitted projections. CBD covers ages 55–90; other methods cover ages 0–100. Population coverage varies by method. Displayed values are mortality rates multiplied by 100. ",a("Methods and research references ↗",href="https://github.com/HG-KaiBai/mortselect#references",target="_blank"))
)
server <- function(input,output,session) {
  selected_data <- reactive({req(input$method);forecasts[[input$method]]})
  output$population_control <- renderUI({
    req(input$mode)
    d <- selected_data(); populations <- sort(unique(d$series$population)); multiple <- input$mode != "population"
    selectizeInput("populations","02  Population",choices=setNames(populations,pretty_population(populations)),selected=if("Canada" %in% populations) "Canada" else populations[1],multiple=multiple)
  })
  output$sex_control <- renderUI({
    req(input$mode)
    if (input$mode == "population") radioButtons("sex","03  Sex",c("Female"="female","Male"="male"),selected="female",inline=TRUE)
    else checkboxGroupInput("sex","03  Sex",c("Female"="female","Male"="male"),selected=c("female","male"),inline=TRUE)
  })
  output$dimension_control <- renderUI({
    req(input$mode)
    d <- selected_data()
    switch(input$mode,
      age=sliderInput("age","04  Age",min=min(d$ages),max=max(d$ages),value=60,step=1),
      population=selectizeInput("ages","04  Ages to compare",choices=d$ages,selected=intersect(c(40,60,80),d$ages),multiple=TRUE),
      year=sliderInput("year","04  Forecast year",min=min(d$years),max=max(d$years),value=max(d$years),step=1,sep=""))
  })
  output$coverage <- renderText({d<-selected_data();sprintf("%d populations · ages %d–%d · %d–%d",length(unique(d$series$population)),min(d$ages),max(d$ages),min(d$years),max(d$years))})
  rows <- reactive({
    d <- selected_data(); shiny::validate(need(length(input$populations)>0,"Choose at least one population."),need(length(input$sex)>0,"Choose at least one sex."))
    ages <- if(input$mode=="year") d$ages else if(input$mode=="population") as.integer(input$ages) else as.integer(input$age)
    years <- if(input$mode=="year") as.integer(input$year) else d$years
    shiny::validate(need(length(ages)>0,"Choose at least one age."),need(length(years)>0,"Choose a forecast year."),need(all(ages %in% d$ages),"Choose an age available for this method."))
    result <- list()
    for (population in input$populations) for (sex in input$sex) {
      i <- which(d$series$population==population & d$series$sex==sex)
      if(!length(i)) next
      matrix <- d$series$values[[i]]
      frame <- expand.grid(age=ages,year=years)
      frame$mortality_rate <- matrix[cbind(match(frame$age,d$ages),match(frame$year,d$years))]
      frame$population <- population;frame$sex <- sex;frame$method <- input$method
      frame$rate_times_100 <- 100*frame$mortality_rate
      frame$series <- if(input$mode=="population") paste("Age",frame$age) else paste(pretty_population(population),tools::toTitleCase(sex),sep=" · ")
      result[[length(result)+1L]] <- frame
    }
    shiny::validate(need(length(result)>0,"No forecasts are available for this selection."))
    do.call(rbind,result)
  })
  output$title <- renderText(switch(input$mode,age=paste("Mortality over time at age",input$age),population="How forecasts vary with age",year=paste("Mortality across ages in",input$year)))
  output$subtitle <- renderText(paste(names(methods)[match(input$method,methods)],"· Saved research forecasts · 2011–2025"))
  chart <- reactive({
    frame <- rows(); frame$x <- if(input$mode=="year") frame$age else frame$year
    p <- ggplot(frame,aes(x=x,y=rate_times_100,color=series,group=series,linetype=sex))+geom_line(linewidth=.85)+
      scale_color_manual(values=rep(c("#247466","#be714b","#6579a3","#a18b3b","#8d6090","#338c9d"),length.out=length(unique(frame$series))))+
      scale_linetype_manual(values=c(female="solid",male="dashed"),guide="none")+
      labs(x=if(input$mode=="year") "Age (years)" else "Forecast year",y="Mortality rate × 100",color=NULL)+
      theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(color="#e9eee7"),legend.position="bottom",axis.text=element_text(color="#697774"),text=element_text(color="#183c38"),legend.text=element_text(size=10))
    if(isTRUE(input$log)){shiny::validate(need(all(frame$rate_times_100>0),"Log scale requires positive rates."));p<-p+scale_y_log10()} else p<-p+scale_y_continuous(limits=c(0,NA))
    p
  })
  output$forecast_plot <- renderPlot(chart(),res=120)
  output$values <- renderTable(rows()[,c("population","sex","age","year","rate_times_100")],digits=6)
  output$download_csv <- downloadHandler(filename=function()paste0("mortselect-",input$method,".csv"),content=function(file)write.csv(rows()[,c("method","population","sex","age","year","mortality_rate","rate_times_100")],file,row.names=FALSE))
  output$download_plot <- downloadHandler(filename=function()paste0("mortselect-",input$method,".png"),content=function(file)ggsave(file,chart(),width=11,height=6,dpi=180,bg="white"))
}
shinyApp(ui,server)
