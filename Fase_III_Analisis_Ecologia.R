# FASE III - DIVERSIDAD FUNCIONAL
# Paramo del Nevado del Ruiz vs. La Tatacoa

# 1. IMPORTAR LA MATRIZ DE RASGOS
datos <- read.csv(
  "matriz_ecologia.csv.csv",
  header = TRUE,
  sep = ",",
  dec = ",",
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# 2. VERIFICAR LOS DATOS
dim(datos)
str(datos)
colSums(is.na(datos))
table(datos$Ecosistema)

# 3. ESTANDARIZAR LOS RASGOS
rasgos <- datos[, 3:6]
rasgos_est <- scale(rasgos)
rownames(rasgos_est) <- datos$Especie

# Calcular Rao por ecosistema



# Calcular distancias funcionales
distancias <- as.matrix(dist(rasgos_est))
distancias <- distancias / max(distancias)

# Definir la función de Rao
calcular_rao <- function(a) {
  p <- a / sum(a)
  as.numeric(t(p) %*% distancias %*% p)
}
rao <- apply(abundancia, 1, calcular_rao)

round(rao, 4)

# Guardar resultados
tabla_rao <- data.frame(
  Ecosistema = names(rao),
  Rao_Q = round(as.numeric(rao), 4)
)

write.csv(tabla_rao,
          "Resultados_Rao.csv",
          row.names = FALSE)

# 6. GRUPOS FUNCIONALES

# Fijar semilla para reproducir los resultados
set.seed(123)

# Agrupar especies según sus cuatro rasgos
grupos <- kmeans(rasgos_est, centers = 3, nstart = 50)

# Asignar un grupo funcional a cada especie
datos$Grupo_funcional <- factor(grupos$cluster)

# Contar especies por grupo y ecosistema
tabla_grupos <- table(datos$Ecosistema,
                      datos$Grupo_funcional)

tabla_grupos


# 7. SHANNON DE GRUPOS FUNCIONALES

# Calcular proporciones de cada grupo
proporciones <- prop.table(tabla_grupos, margin = 1)

# Calcular Shannon
shannon <- -rowSums(
  ifelse(proporciones > 0,
         proporciones * log(proporciones),
         0)
)

round(shannon, 4)

# Guardar resultados
tabla_shannon <- data.frame(
  Ecosistema = names(shannon),
  Shannon_H = round(as.numeric(shannon), 4)
)

write.csv(tabla_shannon,
          "Resultados_Shannon.csv",
          row.names = FALSE)

# 8. ANÁLISIS DE COMPONENTES PRINCIPALES (PCA)

# Ejecutar el PCA con los rasgos estandarizados
pca <- prcomp(rasgos_est,
              center = FALSE,
              scale. = FALSE)

# Ver el resumen del análisis
summary(pca)

# Porcentaje de variación explicado por PC1 y PC2
round(100 * summary(pca)$importance[2, 1:2], 2)


colores <- ifelse(datos$Ecosistema == "Páramo",
                  "darkgreen", "darkorange")

plot(pca$x[,1], pca$x[,2],
     col = colores,
     pch = 19,
     xlab = "PC1 (53.42%)",
     ylab = "PC2 (25.85%)",
     main = "PCA de rasgos funcionales")

abline(h = 0, v = 0, lty = 3, col = "gray70")

legend("topright",
       legend = c("Páramo", "Tatacoa"),
       col = c("darkgreen", "darkorange"),
       pch = 19)


# PASO 13: contribución de los rasgos al PCA

round(pca$rotation[, 1:2], 3)

biplot(pca, scale = 0, cex = 0.65)

# 10. BIPLOT DEL PCA CON GGPLOT2

library(ggplot2)
library(ggrepel)

# Coordenadas de las 20 especies
puntos <- data.frame(
  Especie = datos$Especie,
  Ecosistema = datos$Ecosistema,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

# Contribución de los cuatro rasgos
flechas <- data.frame(
  PC1 = pca$rotation[, 1],
  PC2 = pca$rotation[, 2],
  Rasgo = c("SLA", "LDMC", "Altura", "Densidad")
)

# Ajustar el tamaño visual de las flechas
factor_flechas <- min(
  diff(range(puntos$PC1)) / diff(range(flechas$PC1)),
  diff(range(puntos$PC2)) / diff(range(flechas$PC2))
) * 0.30

flechas$PC1 <- flechas$PC1 * factor_flechas
flechas$PC2 <- flechas$PC2 * factor_flechas

# Porcentajes de variación
varianza <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# Crear la gráfica
grafica_pca <- ggplot(puntos, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray75") +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray75") +
  geom_point(aes(color = Ecosistema), size = 3) +
  geom_text_repel(
    aes(label = Especie, color = Ecosistema),
    size = 2.8,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_segment(
    data = flechas,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.2, "cm")),
    color = "black",
    linewidth = 0.7
  ) +
  geom_text_repel(
    data = flechas,
    aes(x = PC1, y = PC2, label = Rasgo),
    inherit.aes = FALSE,
    fontface = "bold",
    color = "black",
    size = 4
  ) +
  scale_color_manual(
    values = c("Páramo" = "#25854E",
               "Tatacoa" = "#E58A24")
  ) +
  labs(
    title = "Diversidad funcional: Páramo vs. Tatacoa",
    x = paste0("PC1 (", round(varianza[1], 2), "%)"),
    y = paste0("PC2 (", round(varianza[2], 2), "%)"),
    color = "Ecosistema"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Mostrar la gráfica
print(grafica_pca)


# 11. BIPLOT FINAL - VERSION LIMPIA

# Seleccionar especies para destacar
destacadas <- puntos[
  puntos$Especie %in% c(
    "Espeletia hartwegiana",
    "Calamagrostis effusa",
    "Myrsine dependens",
    "Prosopis juliflora",
    "Melocactus curvispinus",
    "Astronium graveolens"
  ), ]

# Crear flechas directamente desde el PCA
vectores <- data.frame(
  Rasgo = c("SLA", "LDMC", "Altura", "Densidad"),
  PC1 = pca$rotation[,1] * 1.8,
  PC2 = pca$rotation[,2] * 1.8
)

# Gráfica
pca_final <- ggplot(puntos, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray80") +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray80") +
  geom_point(aes(color = Ecosistema),
             size = 3.5, alpha = 0.9) +
  geom_text_repel(
    data = destacadas,
    aes(label = Especie, color = Ecosistema),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.6,
    show.legend = FALSE
  ) +
  geom_segment(
    data = vectores,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    inherit.aes = FALSE,
    arrow = grid::arrow(
      length = grid::unit(0.22, "cm")
    ),
    color = "#343A40",
    linewidth = 0.8
  ) +
  geom_label_repel(
    data = vectores,
    aes(x = PC1, y = PC2, label = Rasgo),
    inherit.aes = FALSE,
    fill = "white",
    fontface = "bold",
    size = 3.5,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "Páramo" = "#238B45",
    "Tatacoa" = "#E68A26"
  )) +
  labs(
    title = "Espacio funcional de las comunidades vegetales",
    subtitle = "Páramo del Nevado del Ruiz y La Tatacoa",
    x = "PC1 (53,42 %)",
    y = "PC2 (25,85 %)",
    color = "Ecosistema"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(pca_final)


# 13. PCA - VERSION PARA PUBLICACION

library(ggplot2)
library(ggrepel)

# Porcentaje de variación
var_pca <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# Flechas de los rasgos
vectores <- data.frame(
  Rasgo = c("SLA", "LDMC", "Altura", "Densidad"),
  PC1 = pca$rotation[, 1] * 2.1,
  PC2 = pca$rotation[, 2] * 2.1
)

# Especies que vamos a destacar
etiquetas <- puntos[
  puntos$Especie %in% c(
    "Espeletia hartwegiana",
    "Calamagrostis effusa",
    "Myrsine dependens",
    "Prosopis juliflora",
    "Stenocereus griseus",
    "Astronium graveolens"
  ), ]

# Crear la figura
pca_elegante <- ggplot(puntos, aes(PC1, PC2)) +
  
  # Elipses de dispersión
  stat_ellipse(
    aes(fill = Ecosistema),
    geom = "polygon",
    level = 0.68,
    alpha = 0.10,
    color = NA
  ) +
  
  stat_ellipse(
    aes(color = Ecosistema),
    level = 0.68,
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  
  # Líneas del origen
  geom_hline(yintercept = 0, color = "gray80") +
  geom_vline(xintercept = 0, color = "gray80") +
  
  # Las veinte especies
  geom_point(
    aes(fill = Ecosistema),
    shape = 21,
    size = 4,
    color = "white",
    stroke = 0.7
  ) +
  
  # Etiquetas de especies seleccionadas
  geom_text_repel(
    data = etiquetas,
    aes(label = Especie, color = Ecosistema),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.7,
    point.padding = 0.4,
    max.overlaps = Inf,
    seed = 123,
    show.legend = FALSE
  ) +
  
  # Flechas de los rasgos
  geom_segment(
    data = vectores,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    inherit.aes = FALSE,
    linewidth = 0.9,
    color = "#334155",
    arrow = grid::arrow(
      length = grid::unit(0.25, "cm")
    )
  ) +
  
  # Nombres de los rasgos
  geom_label_repel(
    data = vectores,
    aes(x = PC1, y = PC2, label = Rasgo),
    inherit.aes = FALSE,
    fill = "white",
    color = "#334155",
    fontface = "bold",
    size = 4,
    box.padding = 0.5,
    seed = 123
  ) +
  
  # Paleta de colores
  scale_color_manual(
    values = c("Páramo" = "#20745C",
               "Tatacoa" = "#D98235")
  ) +
  
  scale_fill_manual(
    values = c("Páramo" = "#20745C",
               "Tatacoa" = "#D98235")
  ) +
  
  # Títulos
  labs(
    title = "Espacio funcional de la vegetación",
    subtitle = "Páramo del Nevado del Ruiz y zona semiárida de La Tatacoa",
    x = paste0("PC1 (", round(var_pca[1], 2), " %)"),
    y = paste0("PC2 (", round(var_pca[2], 2), " %)"),
    color = "Ecosistema",
    fill = "Ecosistema",
    caption = "PCA de cuatro rasgos estandarizados | Elipses descriptivas al 68 %"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold", size = 17
    ),
    plot.subtitle = element_text(
      size = 11, color = "gray40"
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      color = "gray50", size = 9
    )
  )

print(pca_elegante)



# 13. PCA - VERSION MÁS LIMPIA

library(ggplot2)
library(ggrepel)

# Porcentaje de varianza explicada
var_pca <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# Tabla de puntos
puntos <- data.frame(
  Especie = datos$Especie,
  Ecosistema = datos$Ecosistema,
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

# Flechas de rasgos
vectores <- data.frame(
  Rasgo = c("SLA", "LDMC", "Altura", "Densidad"),
  PC1 = pca$rotation[,1] * 1.8,
  PC2 = pca$rotation[,2] * 1.8
)

# Solo algunas especies etiquetadas
# Etiquetar las 20 especies
etiquetas <- puntos

pca_elegante <- ggplot(puntos, aes(x = PC1, y = PC2)) +
  stat_ellipse(
    aes(fill = Ecosistema),
    geom = "polygon",
    level = 0.68,
    alpha = 0.10,
    color = NA
  ) +
  stat_ellipse(
    aes(color = Ecosistema),
    level = 0.68,
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_hline(yintercept = 0, color = "gray80") +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_point(
    aes(fill = Ecosistema),
    shape = 21,
    size = 3.5,
    color = "white",
    stroke = 0.6
  ) +
  geom_text_repel(
    data = etiquetas,
    aes(label = Especie, color = Ecosistema),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_segment(
    data = vectores,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black",
    arrow = grid::arrow(length = grid::unit(0.22, "cm"))
  ) +
  geom_label_repel(
    data = vectores,
    aes(x = PC1, y = PC2, label = Rasgo),
    inherit.aes = FALSE,
    fill = "white",
    color = "black",
    fontface = "bold",
    size = 3.5,
    box.padding = 0.3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("Páramo" = "#238B45", "Tatacoa" = "#E67E22")) +
  scale_fill_manual(values = c("Páramo" = "#238B45", "Tatacoa" = "#E67E22")) +
  labs(
    title = "Espacio funcional de la vegetación",
    subtitle = "Páramo del Nevado del Ruiz y zona semiárida de La Tatacoa",
    x = paste0("PC1 (", round(var_pca[1], 2), " %)"),
    y = paste0("PC2 (", round(var_pca[2], 2), " %)"),
    color = "Ecosistema",
    fill = "Ecosistema"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(pca_elegante)


# Guardar PCA con las 20 especies
ggsave(
  "PCA_20_especies.png",
  plot = pca_elegante,
  width = 16,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)


ggsave("PCA_20_especies.png",
       plot = pca_elegante,
       width = 16,
       height = 10,
       dpi = 300)


# 14. GRAFICA DE DIVERSIDAD FUNCIONAL - RAO

grafica_rao <- ggplot(
  tabla_rao,
  aes(x = Ecosistema, y = Rao_Q, fill = Ecosistema)
) +
  geom_col(width = 0.55) +
  geom_text(
    aes(label = sprintf("%.4f", Rao_Q)),
    vjust = -0.6,
    size = 5
  ) +
  scale_fill_manual(values = c(
    "Páramo" = "#238B45",
    "Tatacoa" = "#E67E22"
  )) +
  scale_y_continuous(
    limits = c(0, 0.55)
  ) +
  labs(
    title = "Diversidad funcional de Rao",
    subtitle = "Comparación entre páramo y Tatacoa",
    x = "Ecosistema",
    y = "Entropía cuadrática de Rao (Q)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

print(grafica_rao)


# Guardar figura de Rao

ggsave(
  "Figura_Rao.png",
  plot = grafica_rao,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

file.exists("Figura_Rao.png")


# 15. GRAFICA DE SHANNON FUNCIONAL

grafica_shannon <- ggplot(
  tabla_shannon,
  aes(x = Ecosistema, y = Shannon_H, fill = Ecosistema)
) +
  geom_col(width = 0.55) +
  geom_text(
    aes(label = sprintf("%.4f", Shannon_H)),
    vjust = -0.6,
    size = 5
  ) +
  scale_fill_manual(values = c(
    "Páramo" = "#238B45",
    "Tatacoa" = "#E67E22"
  )) +
  scale_y_continuous(
    limits = c(0, 1.2)
  ) +
  labs(
    title = "Diversidad de grupos funcionales",
    subtitle = "Índice de Shannon (H')",
    x = "Ecosistema",
    y = "Shannon funcional (H')"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

print(grafica_shannon)


# 16. GUARDAR FIGURA DE SHANNON

ggsave(
  "Figura_Shannon.png",
  plot = grafica_shannon,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

file.exists("Figura_Shannon.png")





