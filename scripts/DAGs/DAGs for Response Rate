library(dagitty)
library(ggdag)
library(stringr)
library(dplyr)
library(ggraph)

dag <- dagitty('dag {
  Dose [exposure]
  Response [outcome]
  
  Location -> Stage
  Location -> Dose
  Location -> Response
  
  Age -> Dose
  Age -> Response
  Age -> "Performance Status"
  
  Sex -> Response
  
  "Immune Status" -> Stage
  "Immune Status" -> Agent
  "Immune Status" -> Dose
  "Immune Status" -> Response
  "Immune Status" -> "Performance Status"
  
  Stage -> Agent
  Stage -> Dose
  Stage -> Response
  Stage -> "Performance Status"
  
  Agent -> Dose
  Agent -> Response
  
  "Performance Status" -> Dose
  "Performance Status" -> Response
  
  Year -> Dose
  
  Dose -> Response
}')

dag_tidy <- dag |>
  tidy_dagitty(layout = "circle") |>
  mutate(
    label = str_wrap(name, width = 5),
    edge_type = if_else(name == "Dose" & to == "Response", "main", "other")
  )

dag_dose_response_total <- 
  ggdag(dag_tidy, text = FALSE, node_size = 0) +
  geom_dag_edges(
    aes(edge_colour = edge_type, edge_width = edge_type),
    arrow_directed = arrow(length = unit(10, "pt"), type = "closed")
  ) +
  scale_edge_colour_manual(values = c(main = "red", other = "gray10")) +
  scale_edge_width_manual(values = c(main = 2, other = 0.8)) +
  geom_dag_label(aes(label = label), size = 4) +
  theme_dag() +
  theme(legend.position = "none")

dag_dose_response_total

save_files(
  save_object = dag_dose_response_total,
  filename = "dag_dose_response_total",
  directory = here::here(),
  subD = "files/Directed Acyclic Graphs/Total Effect DAG"
)
# Check adjustment set
adjustmentSets(dag, exposure = "Dose", outcome = "Response")

paths(dag, from = "Dose", to = "Response")

# Analysis 2
## Response ~ Dose + Time + Age + Sex + ImmuneStatus + Location + Stage + Agent + ECOG + Year
#### Time being time to first response assessment
# Analysis 2: Mediation Analysis DAG (Time as mediator)
dag_mediation <- dagitty('dag {
  Dose [exposure]
  Response [outcome]
  Time [mediator]
  
  Location -> Stage
  Location -> Dose
  Location -> Time
  Location -> Response
  
  Age -> Dose
  Age -> Response
  Age -> "Performance Status"
  
  Sex -> Response
  
  "Immune Status" -> Stage
  "Immune Status" -> Agent
  "Immune Status" -> Dose
  "Immune Status" -> Response
  "Immune Status" -> "Performance Status"
  
  Stage -> Agent
  Stage -> Dose
  Stage -> Time
  Stage -> Response
  Stage -> "Performance Status"
  
  Agent -> Dose
  Agent -> Response
  
  "Performance Status" -> Dose
  "Performance Status" -> Time
  "Performance Status" -> Response
  
  Year -> Dose
  Year -> Time
  
  Dose -> Time
  Time -> Response
  Dose -> Response
}')

# Check adjustment set for total effect
adjustmentSets(dag_mediation, exposure = "Dose", outcome = "Response")

# Check adjustment set for direct effect (controlling for mediator)
adjustmentSets(dag_mediation, exposure = "Dose", outcome = "Response", 
               type = "all")

dag_tidy <- dag_mediation |>
  tidy_dagitty(layout = "circle") |>
  mutate(
    label = str_wrap(name, width = 5),
    edge_type = if_else(name == "Dose" & to == "Response", "main", "other")
  )

dag_dose_time_response <- 
  ggdag(dag_tidy, text = FALSE, node_size = 0) +
  geom_dag_edges(
    aes(edge_colour = edge_type, edge_width = edge_type),
    arrow_directed = arrow(length = unit(10, "pt"), type = "closed")
  ) +
  scale_edge_colour_manual(values = c(main = "red", other = "gray10")) +
  scale_edge_width_manual(values = c(main = 2, other = 0.8)) +
  geom_dag_label(aes(label = label), size = 4) +
  theme_dag() +
  theme(legend.position = "none")

dag_dose_time_response


save_files(
  save_object = dag_dose_time_response,
  filename = "dag_dose_time_response",
  directory = here::here(),
  subD = "files/Directed Acyclic Graphs/Mediation Analysis DAG"
)
