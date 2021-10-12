tabel_full <- read.csv('DATASET_def_small.csv', sep=';', header=T)

idx <- which(tabel_full$Country_code_A2=='IS')

add_line <- c(
rowMeans(tabel_full[idx, 5:15], na.rm = T),
rowMeans(tabel_full[idx, 16:27], na.rm = T),
rowMeans(tabel_full[idx, 28:38], na.rm = T),
rowMeans(tabel_full[idx, 39:46], na.rm = T),
rowMeans(tabel_full[idx, 47:57], na.rm = T),
rowMeans(tabel_full[idx, 58:68], na.rm = T),
rowMeans(tabel_full[idx, 69:79], na.rm = T),
rowMeans(tabel_full[idx, 80:90], na.rm = T),
rowMeans(tabel_full[idx, 91:101], na.rm = T),
rowMeans(tabel_full[idx, 102:112], na.rm = T)
)
format(add_line, scientific = FALSE)

colnames(tabel_full)