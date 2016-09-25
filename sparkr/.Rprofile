Sys.setenv(SPARK_HOME = 'spark-1.4.1/')
.libPaths(c(file.path(Sys.getenv('SPARK_HOME'),'R','lib'),.libPaths()))
library(SparkR)
sc <- sparkR.init(master=”local”)