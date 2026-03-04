# custom function retrieved from Stack Overflow for converting ashape to polygon
# https://stackoverflow.com/questions/72449104/converting-alpha-hulls-to-spatial-polygon

hull2poly <- function(my.ashape){
  require(sf)
  if(class(my.ashape) != "ashape") {stop('error, your input must be 
     ashape class')} else
       my.edge<- data.frame(my.ashape$edges)[,c( 'x1', 'y1', 'x2', 'y2')]
     x<- my.edge[,1:2]
     y<- my.edge[,3:4]
     my.edge2<- matrix(t(cbind(x,y)), byrow=T,ncol=2)
     my.edge2<- as.data.frame(my.edge2)
     names(my.edge2)<- c('x','y')
     my.edge2$id <- unlist(lapply((1: (nrow(my.edge2)/2)), 
                                  FUN=function(x){c(rep(x,2))}))
     
     start.edge<- 1
     new.id<- start.edge
     new.edges<- my.edge2[which(my.edge2$id== start.edge ),]
     
     while(length(new.id)<= length(unique(my.edge2$id))-1){
       internal.id<- new.id[length(new.id)]
       edge <- my.edge2[which(my.edge2$id== internal.id ),]
       where.to.search <- my.edge2[which(my.edge2$id %in% new.id ==F ),]
       
       index1<- apply(where.to.search[,1:2], 1, function(x){x == edge[1,1:2]})
       index1<- as.numeric(names(which(apply(index1,2, sum)>0)))[1]
       index2<- apply(where.to.search[,1:2], 1, function(x){x == edge[2,1:2]})
       index2<- as.numeric(names(which(apply(index2,2, sum)>0)))[1]
       main.index<- c(index1, index2)
       
       ifelse(all(!is.na(main.index)), 
              # yes
              {flag<- c(T,T)
              main.index<- main.index[2]
              point.coord<- my.edge2[main.index,] 
              segment<- my.edge2[my.edge2$id==my.edge2[main.index,'id'],]
              new.id<- c( new.id, my.edge2[main.index,]$id) },
              
              # no
              ifelse(which(!is.na(main.index))==1, 
                     # yes
                     {flag<- c(T,F)
                     main.index<- main.index[flag]
                     point.coord<- my.edge2[main.index,] 
                     segment<- 
                       my.edge2[my.edge2$id==my.edge2[main.index,'id'],]
                     new.id<- c( new.id, my.edge2[main.index,]$id)},
                     # no
                     {flag<- c(F,T)
                     main.index<- main.index[flag]
                     point.coord<- my.edge2[main.index,] 
                     segment<- my.edge2[my.edge2$id==my.edge2[main.index,'id'],]
                     new.id<- c( new.id, my.edge2[main.index,]$id)}  ) )
       
       index3<- t(apply(segment, 1, function(x){x ==point.coord}))
       
       new.edges<- rbind(new.edges, rbind(point.coord, segment[which(apply(index3,1, sum)<3),]))
     }
     tst <- st_multipoint(as.matrix(new.edges), dim = "XYZ")
     poly<- tst %>% # 
       st_cast('POLYGON')
     return(poly)}