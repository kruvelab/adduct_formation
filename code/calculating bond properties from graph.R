library(rcdk)
library(DiagrammeR)
library(igraph)
library(rJava)
library(tidyverse)

bond_descriptors <- function(SMILES) {
  #default values
  min_dist_2O = 100
  max_dist_2O = 0
  max_path_non_rotatable_bonds_2O = 0
  min_path_non_rotatable_bonds_2O = 0
  #preprocessign the graph
  m <- parse.smiles(SMILES)[[1]]
  matrix_with_bond_count <- get.connection.matrix(m)
  matrix <- get.adjacency.matrix(m)
  
  atoms <- get.atoms(m)
  names_of_atoms <- tibble(atoms = character())
  for (i in 1:length(atoms)) {
    names_of_atoms <- names_of_atoms %>%
      add_row(atoms = get.symbol(atoms[[i]]))
  }
  names_of_atoms <- names_of_atoms %>%
    mutate(is_O = case_when(
      atoms == "O" ~ TRUE,
      TRUE ~ FALSE),
      is_carbonyl = FALSE
    )
  
  #which of the elements are oxygens
  oxygens <- which(names_of_atoms$is_O) #gives the indexes of O atoms
  for (j in oxygens) {
    #if the oxygen has at least one double bond, it is a carbonyl
    if (max(matrix_with_bond_count[j, ]) > 1) {
      names_of_atoms$is_carbonyl[j] = TRUE
    }
  }
  
  #oxygens----
  #let us measure the number of bonds between oxygens
  #matrix in the graph representation
  graph <- graph_from_adjacency_matrix(matrix)
  #setting the starting conditions for min and max distance
  min_dist_2O = 100
  max_dist_2O = 0
  max_path_non_rotatable_bonds_2O = 0
  min_path_non_rotatable_bonds_2O = 0
  
  if (length(oxygens) > 1) {
    for (i in oxygens) {
      #find all the paths from oxygen i to other elements
      paths <- igraph::shortest_paths(graph, i)
      #let us look specifically at distance to other oxygens
      for (j in oxygens) {
        if (i < j) {
          #path tlenght to the next onygen
          dist <- length(paths$vpath[[j]])
          if (dist > max_dist_2O) {
            max_dist_2O = dist - 1 #-1 because distance to the element itself has been defined as 1 not zero
            max_path_2O = paths$vpath[[j]]
          } 
          if (dist < min_dist_2O) {
            min_dist_2O = dist - 1 #-1 because distance to the element itself has been defined as 1 not zero
            min_path_2O = paths$vpath[[j]]
          }
        }
      }
    }
    #let us calculate the number of non-rotatable bonds between two furthest appart oxygens
    
    path_to_O <- max_path_2O[-1]
    path_to_O <- path_to_O[-length(path_to_O)]
    if (length(path_to_O) > 1) {
      for (k in 1:(length(path_to_O)-1)) {
        first_atom = path_to_O[k]
        next_atom = path_to_O[k+1]
        #let us get the cycles in the graph for this element
        cycles = lapply(all_simple_paths(graph, next_atom, first_atom, mode="out"), function(p) c(first_atom,p))
        cycles = cycles[which(sapply(cycles, length) > 3)]
        if (matrix_with_bond_count[first_atom, next_atom] > 1) {
          max_path_non_rotatable_bonds_2O = max_path_non_rotatable_bonds_2O + 1
        } else if (length(cycles) > 0) {
          max_path_non_rotatable_bonds_2O = max_path_non_rotatable_bonds_2O + 1
        }
      }
    }
    
    #non rotatable bond count for two closest oxygens
    path_to_O <- min_path_2O[-1]
    path_to_O <- path_to_O[-length(path_to_O)]
    if (length(path_to_O) > 1) {
      for (k in 1:(length(path_to_O)-1)) {
        first_atom = path_to_O[k]
        next_atom = path_to_O[k+1]
        #let us get the cycles in the graph for this element
        cycles = lapply(all_simple_paths(graph, next_atom, first_atom, mode="out"), function(p) c(first_atom,p))
        cycles = cycles[which(sapply(cycles, length) > 3)]
        if (matrix_with_bond_count[first_atom, next_atom] > 1) {
          min_path_non_rotatable_bonds_2O = min_path_non_rotatable_bonds_2O + 1
        } else if (length(cycles) > 0) {
          min_path_non_rotatable_bonds_2O = min_path_non_rotatable_bonds_2O + 1
        }
      }
    }
  }
  
  
  
  #carbonyls----
  #now we can also look at the distance to other carbonyls
  carbonyls <- which(names_of_atoms$is_carbonyl) #gives the indexes of =O atoms
  #setting the starting conditions for min and max distance
  min_dist_carbonyl = 100
  max_dist_carbonyl = 0
  max_path_non_rotatable_bonds_2carbonyls = 0
  min_path_non_rotatable_bonds_2carbonyls = 0
  
  if (length(carbonyls) > 1) {
    for (i in carbonyls) {
      #find all the paths from oxygen i to other elements
      paths <- igraph::shortest_paths(graph, i)
      #let us look specifically at distance to other oxygens
      for (j in carbonyls) {
        if (i < j) {
          dist <- length(paths$vpath[[j]])
          if (dist > max_dist_carbonyl) {
            max_dist_carbonyl = dist - 1 #-1 because distance to the element itself has been defined as 1 not zero
            max_path_2carbonyl = paths$vpath[[j]]
          }
          if (dist < min_dist_carbonyl) {
            min_dist_carbonyl = dist - 1 #-1 because distance to the element itself has been defined as 1 not zero
            min_path_2carbonyl = paths$vpath[[j]]
          }
        }
      }
    }
    #let us calculate the number of non-rotatable bonds between two carbonyls
    
    path_to_O <- max_path_2carbonyl[-1]
    path_to_O <- path_to_O[-length(path_to_O)]
    if (length(path_to_O) > 1) {
      for (k in 1:(length(path_to_O)-1)) {
        first_atom = path_to_O[k]
        next_atom = path_to_O[k+1]
        #let us get the cycles in the graph for this element
        cycles = lapply(all_simple_paths(graph, next_atom, first_atom, mode="out"), function(p) c(first_atom,p))
        cycles = cycles[which(sapply(cycles, length) > 3)]
        if (matrix_with_bond_count[first_atom, next_atom] > 1) {
          max_path_non_rotatable_bonds_2carbonyls = max_path_non_rotatable_bonds_2carbonyls + 1
        } else if (length(cycles) > 0) {
          max_path_non_rotatable_bonds_2carbonyls = max_path_non_rotatable_bonds_2carbonyls + 1
        }
      }
    }
    
    path_to_O <- min_path_2carbonyl[-1]
    path_to_O <- path_to_O[-length(path_to_O)]
    if (length(path_to_O) > 1) {
      for (k in 1:(length(path_to_O)-1)) {
        first_atom = path_to_O[k]
        next_atom = path_to_O[k+1]
        #let us get the cycles in the graph for this element
        cycles = lapply(all_simple_paths(graph, next_atom, first_atom, mode="out"), function(p) c(first_atom,p))
        cycles = cycles[which(sapply(cycles, length) > 3)]
        if (matrix_with_bond_count[first_atom, next_atom] > 1) {
          min_path_non_rotatable_bonds_2carbonyls = min_path_non_rotatable_bonds_2carbonyls + 1
        } else if (length(cycles) > 0) {
          min_path_non_rotatable_bonds_2carbonyls = min_path_non_rotatable_bonds_2carbonyls + 1
        }
      }
    }
  }
  
  descriptors <- list("min_dist_2O" = min_dist_2O,
                      "max_dist_2O" = max_dist_2O,
                      "min_dist_carbonyl" = min_dist_carbonyl,
                      "max_dist_carbonyl" = max_dist_carbonyl,
                      "max_path_non_rotatabale_bonds_carbonyl" = max_path_non_rotatable_bonds_2carbonyls,
                      "min_path_non_rotatabale_bonds_carbonyl" = min_path_non_rotatable_bonds_2carbonyls,
                      "max_path_non_rotatabale_bonds_2O" = max_path_non_rotatable_bonds_2O,
                      "min_path_non_rotatabale_bonds_2O" = min_path_non_rotatable_bonds_2O)
  return(descriptors)
}


test <- bond_descriptors("COC1C(CC2CN3CCC4=C(C3CC2C1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC")
