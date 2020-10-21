library(here)
library(readxl)
library(dplyr)
library(igraph) # basic package for network analysis in R
library(visNetwork) # interactive visualisations
library(xlsx)
library(ggplot2)
library(stargazer)
library(sjPlot)

#---Read in and reformat publication list----

soc_pub <- as.data.frame(read_excel("soc_pub_all-sci.xlsx"))
head(soc_pub, 6)

soc_pub.coauth <- sapply(as.character(soc_pub$author), strsplit, "; ") %>% 
  lapply(trimws)

soc_pub.coauth.unique <- unique(unlist(soc_pub.coauth))[order(unique(unlist(soc_pub.coauth)))]

#---Creating Edge Set----

soc_pub.bipartite.edges <- lapply(soc_pub.coauth, function(x) {soc_pub.coauth.unique %in% x})

soc_pub.bipartite.edges <- do.call("cbind", soc_pub.bipartite.edges) # dimension is number of authors x number of papers

rownames(soc_pub.bipartite.edges) <- soc_pub.coauth.unique

mat <- soc_pub.bipartite.edges %*% t(soc_pub.bipartite.edges)

soc_pub.mat <- mat[order(rownames(mat)), order(rownames(mat))]

#---Create Network Object----

soc_pub.net <- graph_from_adjacency_matrix(as.matrix(soc_pub.mat), mode = "undirected", weighted = T, diag = F) 

#---Read Author Attributes----

auth_att <- as.data.frame(read_excel("auth_attributes.xlsx"))
auth_att <- auth_att %>%
  mutate(dep_cat = if_else(is.na(dep_cat), 6, dep_cat), # 6 = external author
         socium = if_else(dep_cat == 6, 2, 1)) 

# merge attribute data to graph
df.soc_pub <- igraph::as_data_frame(soc_pub.net, 'both')  # decompose graph into nodelist and edgelist
df.soc_pub$vertices <- df.soc_pub$vertices %>%            # add node attributes
  left_join(auth_att) %>%
  mutate(id = c(1:length(name)))

soc_pub.net.att <- graph_from_data_frame(df.soc_pub$edges, directed = F, vertices = df.soc_pub$vertices)

#---Computing Network Measures----

n <- soc_pub.net.att

# Measures for individual nodes
V(n)$degree <- degree(n, normalized = FALSE)
V(n)$betweenness <- betweenness(n, normalized = FALSE)
V(n)$eigenvector <- as.numeric(eigen_centrality(n)$vector)
V(n)$strength <- strength(n)

# Measures for overall network
n$comp <- components(n)
n$ncomp <- count_components(n)
n$centdeg <- centralization.degree(n)
n$centbtw <- centralization.betweenness(n)
n$centeig <- centralization.evcent(n)

# Prep for visualisations
V(n)$size <- 1+0.5*degree(n)
E(n)$width <- E(n)$weight
V(n)$label.cex <- .7
V(n)$color <- c("orange", "yellow", "red", "green", "purple")[V(n)$dep_cat]
E(n)$color <- "gray80"

#---Subnetwork of Socium members----

n_sub <- delete_vertices(n, V(n)[dep_cat == 6]) # delete vertices that are external co-authors

E(n_sub)$width <- 0.8*E(n_sub)$weight
V(n_sub)$size <- 5+0.1*(degree(n_sub))^2
n_sub$assort <- assortativity.nominal(n_sub, type = V(n_sub)$dep_cat) # department homophily

write.graph(n_sub, file = here("visone_graphs", "n_sub.graphml"), format = "graphml") # formatting for visone

#---Interactive Visualisation with visNetwwork----

visIgraph(n_sub)

#---Extracting Network Measures----
#--Macro Descriptives
#-complete network

authors <- vcount(n)
publications <- length(soc_pub$title)
edges <- ecount(n)

pub_per_auth <- mean(rowSums(soc_pub.bipartite.edges))
auth_per_pub <- mean(colSums(soc_pub.bipartite.edges))

#-socium sub-network

authors_soc <- vcount(n_sub)
edges_soc <- ecount(n_sub)

#-put into table
macro_des <- data.frame("1" = c("No. Publications", "No. Authors", "No. Edges", "Publications/Author", "Authors/Publication",
                                 "Degree Centralization", "Betweenness Centralization", "Eigenvector Centralization",
                                 "No. Components", "No. of Authors in Main Component (%)", "Isolates",
                                 "No. Authors", "No. Edges", "Department Assortativity"), # socium
                        "2" = (c(publications, authors, edges, round(c(pub_per_auth, auth_per_pub, 
                                 n$centdeg$centralization, n$centbtw$centralization, n$centeig$centralization), digits = 2), 
                                 n$ncomp, "1081 (87.2%)", sum(degree(n)==0),
                                 authors_soc, edges_soc, round(n_sub$assort, digits = 2))))

save(macro_des, file = "macro_des.RData")
#write.xlsx(macro_des, "macro_des.xlsx", row.names = F, col.names = F)

#--Micro Descriptives

df.vert_n <- igraph::as_data_frame(n, 'vertices')
df.vert_red <- filter(df.vert_n, dep_cat != 6) # reduced network (only SOCIUM members, but with measures from orig. netw.)
save(df.vert_red, file = "micro_des_soc.RData")

stargazer(df.vert_red[c("dep1", "dep2", "dep3", "dep4", "dep5", "dep6", "dep_cat2", "tot_pub", "degree", "betweenness", "eigenvector", "strength")], 
          title = "Summary statistics of author characteristics", digits = 1, type = "text", out = "summary.htm")

(top20_deg <- df.vert_red %>%
    arrange(desc(degree)) %>%
    select(name, degree, dep_cat) %>%
    head(20))

(top20_btw <- df.vert_red %>%
    arrange(desc(betweenness)) %>%
    mutate(betweenness = round(betweenness, 1)) %>%
    select(name, betweenness, dep_cat) %>%
    head(20))

(top20_eig <- df.vert_red %>%
    arrange(desc(eigenvector)) %>%
    mutate(eigenvector = round(eigenvector, 1)) %>%
    select(name, eigenvector, dep_cat) %>%
    head(20))

(top20_str <- df.vert_red %>%
    arrange(desc(strength)) %>%
    mutate(strength = round(strength, 1)) %>%
    select(name, strength, dep_cat) %>%
    head(20))

top20 <- data.frame(top20_deg, top20_btw, top20_eig, top20_str)

save(top20, file = "top20_df.RData")
#write.xlsx(top20, "top20.xlsx")

#---Further Anylyses----

#-H3: Degree of health department

summary(df.vert_red$degree)

ggplot(df.vert_red, aes(x=degree)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks = seq(0, 175, 10))

deg_dat <- df.vert_red %>%
  filter(degree>0)

ggplot(deg_dat, aes(x=degree)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks = seq(0, 175, 10))

deg_dat <- deg_dat %>%
  mutate(lndeg = log(degree)) %>%
  mutate(dep_cat = factor(dep_cat, labels = c("norm", "econ", "ineq", "life", "health")))

deg_dat_rem <- deg_dat %>%
  filter(!(id == 330 | id == 881))

ggplot(deg_dat_rem, aes(x=lndeg)) +
  geom_histogram() 

summary(m1 <- lm(lndeg ~ factor(dep_cat), data = deg_dat_rem))
save(m1, file = "reg_results.RData")

stargazer(m1, star.cutoffs = c(0.05, 0.01, 0.001), type = "text", out = "m1.htm",
          covariate.labels = c("Economy", "Inequality", "Life Course", "Health"),
          dep.var.labels = "ln(Degree)", dep.var.caption = "")

# Coef-Plot

p <- plot_model(m1, type = "pred", terms = "dep_cat", 
                title = "Predicted Values of Degree Centrality", axis.title = c("Department", "Degree"))

p + theme_sjplot2()

