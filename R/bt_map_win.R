bt_map_win <- function(region, record,
                       mode = c("count","sum","mean","all")) {

  mode <- match.arg(mode)

  value_col <- names(record)[ncol(record)]

  record[[value_col]] <- as.numeric(record[[value_col]])

  region <- as.data.frame(region)
  record <- as.data.frame(record)

  names(region)[1:3] <- c("chr","start","end")
  names(record)[1:3] <- c("chr","start","end")

  region$region_id <- seq_len(nrow(region))

  region_dt <- data.table::as.data.table(region)
  record_dt <- data.table::as.data.table(record)

  data.table::setkey(region_dt, chr,start,end)
  data.table::setkey(record_dt, chr,start,end)

  hits <- data.table::foverlaps(
    record_dt,
    region_dt,
    nomatch = 0L
  )

  if(nrow(hits)==0){

    region$overlap_n <- 0

    if(mode %in% c("sum","all"))
      region$overlap_sum <- 0

    if(mode %in% c("mean","all"))
      region$overlap_mean <- 0

    return(region)
  }

  stat_dt <- data.frame(
    region_id = hits$region_id,
    value = as.numeric(hits[[value_col]])
  )

  count_df <- aggregate(
    value ~ region_id,
    stat_dt,
    length
  )

  names(count_df)[2] <- "overlap_n"

  result <- merge(
    region,
    count_df,
    by="region_id",
    all.x=TRUE
  )

  if(mode %in% c("sum","all")){

    sum_df <- aggregate(
      value ~ region_id,
      stat_dt,
      sum,
      na.rm=TRUE
    )

    names(sum_df)[2] <- "overlap_sum"

    result <- merge(
      result,
      sum_df,
      by="region_id",
      all.x=TRUE
    )
  }

  if(mode %in% c("mean","all")){

    mean_df <- aggregate(
      value ~ region_id,
      stat_dt,
      mean,
      na.rm=TRUE
    )

    names(mean_df)[2] <- "overlap_mean"

    result <- merge(
      result,
      mean_df,
      by="region_id",
      all.x=TRUE
    )
  }

  result[is.na(result)] <- 0

  result[order(result$region_id),]
}
