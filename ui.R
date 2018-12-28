
fluidPage(
	titlePanel("TCGA Gene expression & survival"),

# p("A webapp to interrogate the curated TCGA clinical data."),


	sidebarPanel(
		textInput(inputId="gene",
			value="MET",
			label="Enter gene symbol"
			),
		radioButtons("type", "Survival type:",
			c("Overall Survival (OS)" = "OS",
				"Progression Free Survival" = "PFS")),

		conditionalPanel(condition="input.tabselect==1",
			sliderInput(inputId="hvslowcut",label="High versus low cut",value=50,min=10,max=90,step=5)
			),

		conditionalPanel(condition="input.tabselect==2",
			sliderInput(inputId="intervalcut",label="Lower and Upper cut",value=c(25,75),min=10,max=90,step=5)
			),
		br(),
		# fluidRow(
		# 	column(4,htmlOutput("HRtxt")),
		# 	column(1, plotOutput("HR") )
		# 	)
		h5("Survival association to the continuous expression value:"),
		tableOutput("HRtxt"),
		downloadButton("downloadData", "Download gene data")

		),

	mainPanel(
		tabsetPanel(
			tabPanel("Split", value=1,plotOutput("hvslowplot")), 
			tabPanel("Interval", value=2,plotOutput("intervalplot")), 
			# tabPanel("Data table", value=3,tableOutput("table")),
			tabPanel("Info", value=3,
				p("This web application proposes an easy to use interface to interogate a curated version of the TCGA pancreatic cancer dataset (PAAD)."),
				p("The data shown in this app only uses the verified pancreatic ductal adenocarcinoma after the exclusion of pancreatitis, neuroendocrine tumors and other irrelevant samples. This sums to 150 samples for Overall survival analysis (OS) and 144 for Progression Free Survival analysis (PFS)."),
				p("The clinical and transcriptomic data were taken from the TCGA repository as described in:"),
				p("Raphael, B. J., Hruban, R. H., Aguirre, A. J., Moffitt, R. A., Yeh, J. J., Stewart, C., ... & Gabriel, S. B. (2017). Integrated genomic characterization of pancreatic ductal adenocarcinoma. Cancer Cell, 32(2), 185-203."),
				br(),

				p("This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): http://cit.liguecancer.net. For any question please contact CITR@ligue-cancer.net")
				),
			id="tabselect"
			)
		

		)

	)
