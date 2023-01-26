

# x = buffer_pol_union
# max_dist = max_distance_between_pts
# near_buffer_dist = near_buffer_distance

#' Estimates the midline(s) of sf polygon(s)
#'
#' Uses Voronoi polygons to estimate the midlines of one or more sf polygons
#'
#' Taking an sf polygon or feature collection of polygons, the function uses Voronoi tessellation to estimate the polygon midlines. Sufficient density of points are required on the polygon boarder facilitate the Voronoi tessellation. Large gaps between points are possible on straight lines on polygon perimeters, therefore dfMaxLength is used to stipulate the maximum distance between points and add points where required. This argument is passed to \code{\link[sf:geos_unary]{sf::st_segmentize()}}.
#'
#' The Voronoi tessellation is likely to lead to extraneous lines which do not form part of the intended midline(s). Additional functions \code{\link{midlines_clean}} and \code{\link{midlines_check}} will hopefully help to deal with these.
#'
#' Where there is a region of interest defined by an sf linestring, e.g. of a bounding box, this can be specified to ensure the midlines do not extend beyond this.
#'
#' @param x an sf polygon within which to estimate the midline
#' @param border_line an sf linestring forming the exterior border of the area of interest
#' @param dfMaxLength maximum distance between points in polygon x used to generate Voronoi polygons. Argument passed to \code{\link[sf:geos_unary]{sf::st_segmentize()}}.
#'
#' @examples
#' library(sf)
#' poly = st_buffer(st_linestring(matrix(c(0,0,10,0,10,10,0,10,0,0),ncol=2, byrow=TRUE) ),0.75)
#' plot(poly, col = "GRAY")
#'
#' ml = midlines_draw(poly, dfMaxLength = 1)
#' plot(ml$geometry, add = TRUE)
#'
#' @export
midlines_draw <- function(polygon, border = NULL, max_length = NULL){
  
  if(!(any(class(st_geometry(polygon)) == "sfc_POLYGON"))){
    stop("Input must be a valid sfc polygon.")
  }

  if (!is.null(max_length) && !is.numeric(max_length)){
    stop("max_length must be a numeric value.")
  }

  if (!is.null(border)){
    stop("border must be a valid sfc object.")
  }
 
  multiline = st_union(st_cast(polygon, "MULTILINESTRING"))
 
  if (!is.null(max_length)){
    multiline = st_segmentize(multiline, dfMaxLength = max_length)
  }
 
  voronoi_edges = do.call("c", st_geometry(multiline))  |>
    st_voronoi(bOnlyEdges = TRUE) |>
    st_sfc() |>
    st_sf(crs = sf::st_crs(polygon)) |>
    st_cast("LINESTRING")
  st_geometry(voronoi_edges) <- "geometry"
  
  voronoi_edges = voronoi_edges[unlist(st_contains_properly(polygon, voronoi_edges)),]
  voronoi_edges$line_id = 1:nrow(voronoi_edges)

  if (!is.null(border)) {
    border_poly = st_sfc(st_polygon(st_union(border)), crs = st_crs(border))
    voronoi_edges = st_intersection(voronoi_edges, border_poly)
  }
  
  return(voronoi_edges)
}






# x = midlines_all
# n_removed=10
# border_line = bbox_as_line
# border_distance = set_units(1,"m")
# i=1

#' Aims to clean extraneous lines from estimated midlines
#'
#' Intended for use following \code{\link{midlines_draw}} which uses Voronoi tessellation to estimate polygon midlines. The Voronoi tessellation results in extraneous lines, in addition to the intended midlines. This function aims to remove those lines.
#'
#' Extraneous lines are often short deadends protruding from the intended midlines. This function identifies these lines by identifying line ends and flagging them (with the addition of a 'flag' variable). The process iterates through several cycles of line end identification with the number of cycles specified by the option n_removed. It is likely that some of the intended midlines will also be flagged, it their ends. All lines are returned so that the user can clearly see which lines have been flagged and all lines can then be passed to \code{\link{midlines_check}} which should enable the midlines which have been flagged to be differentiated from the extraneous lines. Depending on the specific use, it may be best to use this function (\code{\link{midlines_check}}) more than once.
#'
#' The border_lines option preemptively prevents lines being flagged if they intersect with a boarder defined by an sf linestring. This might be useful if the border intersects with the extremities of the midlines to prevent their being flagged for removal.
#'
#' @param x Simple features collection. Intended for use with the output from \code{\link{midlines_draw}}
#'
#' @param n_removed specified the number of cycles of line removal
#'
#' @param border_line an sf linestring forming the exterior border of the area of interest (see below)
#'
#' @examples
#' library(sf)
#' # 1
#' poly = st_buffer(st_linestring(matrix(c(0,0,10,0,10,10,0,10,0,0),ncol=2, byrow=TRUE) ),0.75)
#' plot(poly, col = "GRAY")
#'
#' ml = midlines_draw(poly, dfMaxLength = 1)
#' plot(ml$geometry, add = TRUE)
#'
#' ml_clean = midlines_clean(ml)
#' plot(ml_clean$geometry, col = ml_clean$removed_flag, add = TRUE)
#'
#' #2
#' p1 = st_buffer(st_linestring(matrix(c(0,0,30,0),ncol=2, byrow=TRUE) ),0.75)
#' plot(p1)
#' p2 = st_buffer(st_linestring(matrix(c(9,5,9,0,20,0,18,-4),ncol=2, byrow=TRUE) ),0.75)
#' plot(p2, add = TRUE)
#' p3 = st_union(p1, p2)
#' plot(p3, col = "GRAY")
#'
#' bbox_as_line = st_cast(st_as_sfc(st_bbox
#'   (c(xmin = 0, xmax = 30, ymax = -10, ymin = 10))),"LINESTRING")
#' plot(bbox_as_line, add = TRUE)
#'
#' ml = midlines_draw(p3, dfMaxLength = 1)
#' plot(ml$geometry, add = TRUE)
#'
#' ml_clean = midlines_clean(ml, n_removed = 10)
#' plot(ml_clean$geometry, col = ml_clean$removed_flag, add = TRUE)
#'
#' ml_clean2 = midlines_clean(ml, n_removed = 10, border_line = bbox_as_line)
#' plot(p3, col = "GRAY")
#' plot(ml_clean2$geometry, col = ml_clean2$removed_flag, add = TRUE)
#'
#' @export
midlines_clean = function(x, n_removed = 1, border_line = NULL){

  # Identify those edges that intersect with the borderline (if specified)
  if(!(is.null(border_line))) {
    x$border_intersect = as.vector(sf::st_intersects(x$geometry, border_line , sparse = FALSE))
  }

  if (n_removed < 1) {
      stop("n_removed must be a positive integer greater than zero.")
    }
  # to prevent the warning about repeat attributes for all sub-geometries
  st_agr(x) = "constant"

  mid_points = sf::st_cast(x, "POINT")


  mid_points$point_id = 1:nrow(mid_points)

  for(i in 1:n_removed) {
    # Identify those edges that are of length 1
    mid_points |> 
      group_by(geometry) |>
      mutate(dead_point = n() <= 1) -> mid_points

    if(!is.null(border_line)) {
      mid_points$dead_point[mid_points$border_intersect==TRUE] <- FALSE
    }

    mid_points <- mid_points |>
      group_by(line_id) |>
      filter(!any(dead_point))
  }

  trimmed_lines = mid_points %>%
      group_by(line_id) %>%
      summarise(do_union = FALSE) %>%
      st_cast("LINESTRING")

  return(trimmed_lines)
}

group_id = line_id = geometry = NULL


# plot(thames, add = TRUE, lc = "blue"); plot(test |> midlines_clean(n_removed = 4) |> midlines_debit(length = units::set_units(5000, "m")) |> midlines_smooth(width = 10), add = TRUE, col = "red")
#  test[1,]$geometry[[1]][1,] = from 
#  test[1,]$geometry[[1]][2,] = to


trim_ends <- function(graph, n = 1) {
  graph |> activate(nodes) -> graph
  for (i in 1:n){
     graph <- graph |> mutate(deg = centrality_degree()) |> filter(deg > 1)
  }
  return(graph |>select(-deg))
}


midlines_draw(thames) |> 
  #midlines_debit(length = units::set_units(5000, "m")) |>
  mutate(length = st_length(geometry)) -> midline


midline |> group_by(line_id) |>
  mutate(
    from = list(st_point(geometry[[1]][1,])),
    to   = list(st_point(geometry[[1]][2,])))|>
  filter(!is.infinite(length)) |> as_tibble() |> select(-geometry) -> edges


igraph::graph_from_data_frame(edges |> relocate(from, to),  directed = FALSE) -> gr

tidygraph::as_tbl_graph(gr)  |> tidygraph::activate(edges) -> grph

grph |> igraph::distances(weights = as_tibble(grph)$length) -> distance_edges

distance_edges[is.infinite(distance_edges)] <- 0
which(distance_edges == max(distance_edges), arr.ind = TRUE) -> from_to
from_to
shortest_paths(
  grph,
  from = from_to[1],
  to = from_to[2],
  weights = as_tibble(grph)$length,
  output = "vpath")$vpath[[1]] -> shortest




grph |> 
  activate(nodes) |> 
  filter(name %in% names(shortest)) |> 
  trim_ends(5) |>
  mutate(is_end = centrality_degree() == 1) |>
  arrange(factor(name, levels =  as.factor(names(shortest)))) |>
  activate(edges) |>
  mutate(is_end = .N()$is_end[from] | .N()$is_end[to]) |> select(line_id, is_end) |> as_tibble() -> line_to_keep


midline |> 
  inner_join(line_to_keep, by = "line_id") |> arrange(from)  -> medial_axis




extended_lines <- (medial_axis$geometry[which(medial_axis$is_end)] - st_centroid(medial_axis$geometry[which(medial_axis$is_end)])) * 20 + st_centroid(medial_axis$geometry[which(medial_axis$is_end)])

st_crs(extended_lines) <- st_crs(thames)
end_lines <- st_intersection(extended_lines, thames)

point_to_change <- which(rowSums(extended_lines[[1]][,] - end_lines[[1]][,]) != 0)
medial_axis$geometry[which(medial_axis$is_end)][[1]][point_to_change,]  <- end_lines[[1]][point_to_change,]
point_to_change <- which(rowSums(extended_lines[[2]][,] - end_lines[[2]][,]) != 0)
medial_axis$geometry[which(medial_axis$is_end)][[2]][point_to_change,]  <- end_lines[[2]][point_to_change,]





CreateSegments <- function(coords, length = 0, n.parts = 0) {
    stopifnot((length > 0 || n.parts > 0))
    # calculate total length line
    total_length <- 0
    for (i in 1:(nrow(coords) - 1)) {
        d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i + 
            1, 2])^2)
        total_length <- total_length + d
    }

    # calculate stationing of segments
    if (length > 0) {
        stationing <- c(seq(from = 0, to = total_length, by = length), total_length)
    } else {
        stationing <- c(seq(from = 0, to = total_length, length.out = n.parts), 
            total_length)
    }

    # calculate segments and store the in list
    newlines <- list()
    for (i in 1:(length(stationing) - 1)) {
        newlines[[i]] <- CreateSegment(coords, stationing[i], stationing[i + 
            1])
    }
    return(newlines)
}

CreateSegment <- function(coords, from, to) {
    distance <- 0
    coordsOut <- c()
    biggerThanFrom <- F
    for (i in 1:(nrow(coords) - 1)) {
        d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i + 
            1, 2])^2)
        distance <- distance + d
        if (!biggerThanFrom && (distance > from)) {
            w <- 1 - (distance - from)/d
            x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
            y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
            coordsOut <- rbind(coordsOut, c(x, y))
            biggerThanFrom <- T
        }
        if (biggerThanFrom) {
            if (distance > to) {
                w <- 1 - (distance - to)/d
                x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
                y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
                coordsOut <- rbind(coordsOut, c(x, y))
                break
            }
            coordsOut <- rbind(coordsOut, c(coords[i + 1, 1], coords[i + 1, 
                2]))
        }
    }
    return(coordsOut)
}

medial_axis |> midlines_smooth(10) -> smooth_medial_axis

c(tail(st_geometry(smooth_medial_axis),-1), st_geometry(smooth_medial_axis)[1]) |> st_sf() |> mutate(line_id = row_number()) -> medial_axis

st_geometry(medial_axis) |>
  st_cast("POINT") |>
  unique() |>
  st_sfc() |>
  unlist() |>
  matrix(ncol = 2, byrow = TRUE) |>
  as.data.frame() |>
  CreateSegments(n.parts = 100)  -> segments
purrr::map(segments,\(mat){mat[nrow(mat),] |> st_point()}) |> st_sfc() -> split_pts



st_crs(split_pts) <- st_crs(medial_axis)
n <- 20
purrr::map(st_geometry(medial_axis), \(line) {
  purrr::map2(seq(0,1-1/n,1/n), seq(1/n,1,1/n), \(from, to){lwgeom::st_linesubstring(line, from, to)}) |> st_sfc() |> st_sf()
}) -> list_sf
lines <- bind_rows(list_sf)
st_geometry(lines) <- "geometry"





st_crs(lines) <- st_crs(medial_axis)
st_geometry(medial_axis)[st_nearest_feature(medial_axis, split_pts)] 

st_geometry(lines)[st_nearest_feature(split_pts, lines)] 

st_geometry(lines)[st_nearest_feature(split_pts, lines)] |> 
st_snap(split_pts, tolerance = 0) -> snapped 
snapped + (snapped |> st_cast("POINT"))[c(FALSE, TRUE)] - st_centroid(snapped) -> snapped_bis
snapped + (snapped |> st_cast("POINT"))[c(TRUE, FALSE)] - st_centroid(snapped) -> snapped

st_crs(snapped_bis) <- st_crs(medial_axis)
st_crs(snapped) <- st_crs(medial_axis)
# c(
# snapped[which(abs(st_distance(split_pts, snapped)) , arr.ind = TRUE)[,1]],
# snapped_bis[which(abs(st_distance(split_pts, snapped_bis)), arr.ind = TRUE)[,1]]) -> good_snap

apply(
matrix(c(
apply(abs(st_distance(split_pts, snapped)), 2, min),
apply(abs(st_distance(split_pts, snapped_bis)), 2, min)), ncol = 2), 1, which.min) -> pos_min

c(snapped[pos_min == 1], snapped_bis[pos_min == 2])  -> good_snap

good_snap |>
  purrr::map(\(line){
    st_rotate_around(
      line,
      theta = pi/2,
      around = st_rotate_around(st_centroid(line)))
    }) |> st_sfc()  -> perp_lines


ext_perp_lines <- (perp_lines - st_centroid(perp_lines)) * 300 + st_centroid(perp_lines)
st_crs(ext_perp_lines) <- st_crs(medial_axis)



ggplot() +
  geom_sf(data = thames) +
  geom_sf(data = medial_axis) +
  geom_sf(data = split_pts, color = "red") +
  #geom_sf(data = good_snap, color = "blue") + 
  geom_sf(data = ext_perp_lines, color = "green")

 st_split(thames, ext_perp_lines) |> st_collection_extract(c("POLYGON")) -> splited_medial_axis

splited_medial_axis |>
  mutate(id =row_number()) |>
  ggplot() +
  geom_sf(aes(geometry = geometry, fill = id))





st_rotate_around <- function(geometry, theta = 0, around = sf::st_point(c(0, 0))) {
  rotation_matrix <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  rotated_geometry <- ((geometry - around) * rotation_matrix) + around
  # if (sf::st_geometry_type(rotated_geometry) == "POLYGON") {
  #   colnames(rotated_geometry[[1]]) <- c("x", "y")
  # }
  return(rotated_geometry)
}


data(kola.background)
library(mvoutlier)

coast <- kola.background$coast[147:351, ]

class(coast)

odd <- function(x) x%%2 != 0

even <- function(x) x%%2 == 0



st_geometry(medial_axis[[1]])


splitLines<- function(spobj, dist, start = T, sf = F){
    xydf<-coordBuild(spobj)
    if (start == F){
        xydf<-xydf[rev(rownames(xydf)),]
    }
    spoints <- split(xydf, 100)
    linelist <- list()
    lineslist <- list()
    id <- 1
    if(!sf) {
        j <- 1
        for(i in 1:(nrow(spoints)-1)){
            linelist[j] <- Line(spoints[c(i, i + 1), c(1:2)])
            j = j + 1
            if(spoints[i+1,3] == 1){ 
                lineslist[id]<-Lines(linelist, ID = id)
                id = id+1
                linelist<-list()
                j = 1
            }
        }
        return(SpatialLinesDataFrame(SpatialLines(lineslist), data = data.frame(id = 0:(length(lineslist)-1))))
    } else {
        start <- 1
        for(i in 1:(nrow(spoints)-1)){
            if(spoints[i+1,3] == 1){ 
                lineslist[[id]] <- sf::st_linestring(as.matrix(spoints[c(start:(i + 1)), c(1:2)], ncol = 2))
                id <- id + 1
                start <- i + 1
            }
        }
        return(sf::st_sf(id = 1:length(lineslist), geom = sf::st_sfc(lineslist)))
    }
}


coordBuild <- function(spobj){
    if("sfc_LINESTRING" %in% class(spobj)) {
      return(setNames(data.frame(sf::st_zm(sf::st_coordinates(spobj))), c("x", "y")))
    }
    if(class(spobj) %in% c("SpatialLinesDataFrame",    "SpatialLines")){
        coords <- lapply(spobj@lines, function(x) lapply(x@Lines, function(y) y@coords))
        coords <- ldply(coords, data.frame)
        names(coords) <- c("x","y")
        return(coords)
    }
    if(class(spobj) %in% c("SpatialPolygonsDataFrame", "SpatialPolygons")){
        coords <- lapply(spobj@polygons, function(x) lapply(x@Polygons, function(y) y@coords))
        coords <- ldply(coords, data.frame)
        names(coords) <- c("x","y")
        return(coords)
    }
    if(class(spobj) == "data.frame"){
        if(all(c("x", "y") %in% tolower(names(spobj)))){
            return(spobj[c("x", "y")])
        }
        else{
            stop("Dataframe provided does not have x, y columns")
        }
    }
    stop("Class of spatial argument is not supported. Need SpatialLinesDataFrame or SpatialPolygonsDataFrame")
}
