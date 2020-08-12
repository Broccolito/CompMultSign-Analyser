center_mean = function(x,k = 3){
  # k is the kernal size
  xs = vector()
  lx = length(x)
  span = floor(k/2)
  for(i in 1:lx){
    from = ifelse((i - span)>=1,(i-span),1)
    to = ifelse((i + span)<=lx,(i+span),lx)
    kernal = x[from:to]
    xs = c(xs,mean(kernal,na.rm = TRUE))
  }
  return(xs)
}

exp_base = function(x,base = exp(1)){
  return(base^x)
}

process_signal = function(x,exponential = exp(1),
                          kernal_size = 3,
                          smooth_method = "3RS3R"){
  x = x %>% 
    exp_base(base = exponential) %>%
    center_mean(k = kernal_size) %>%
    smooth(kind = smooth_method) %>%
    as.numeric()
  return(x)
}

peakdet = function(v, delta, x = NULL){
  maxtab = NULL
  mintab = NULL
  if(is.null(x)){
    x = seq_along(v)
  }
  if(length(v) != length(x)){
    stop("Input vectors v and x must have the same length")
  }
  if(!is.numeric(delta)){
    stop("Input argument delta must be numeric")
  }
  if(delta <= 0){
    stop("Input argument delta must be positive")
  }
  mn = Inf
  mx = -Inf
  mnpos = NA
  mxpos = NA
  lookformax = TRUE
  for(i in seq_along(v)){
    this = v[i]
    if(this > mx){
      mx = this
      mxpos = x[i]
    }
    if(this < mn){
      mn = this
      mnpos = x[i]
    }
    if(lookformax){
      if (this < mx - delta){
        maxtab = rbind(maxtab, data.frame(pos = mxpos, val = mx))
        mn = this
        mnpos = x[i]
        lookformax = FALSE
      }
    }
    else{
      if(this > mn + delta){
        mintab = rbind(mintab, data.frame(pos = mnpos, val = mn))
        mx = this
        mxpos = x[i]
        lookformax = TRUE
      }
    }
  }
  return(list(maxtab = maxtab, mintab = mintab))
}

find_peak = function(x,threshold_percentage = 0.95){
  peak_df = peakdet(x, delta = quantile(x,threshold_percentage))
  peak_df = peak_df$maxtab
  return(peak_df)
}

generate_processed_df = function(chromo=1,
                                 exponential = exp(1),
                                 kernal_size = 3,
                                 smooth_method = "3RS3R",
                                 threshold = FALSE,
                                 ranking = TRUE,
                                 threshold_percentage = 0.95,
                                 top_ranking = 50){
  cms = filter(cms_list,chr == chromo)$cms_score
  loc = filter(cms_list,chr == chromo)$loc
  processed_signal = process_signal(cms,
                                    exponential = exponential,
                                    kernal_size = kernal_size,
                                    smooth_method = smooth_method)
  processed_df = data.frame(loc = loc, 
                            score = processed_signal,
                            peak = 0)
  if(threshold){
    peaks = find_peak(processed_signal,threshold_percentage)
    processed_df[,3][peaks$pos] = 1
  }
  if(ranking){
    peaks = find_peak(processed_signal,0)
    processed_df[,3][peaks$pos] = 2
    processed_df = arrange(processed_df,desc(score))
    temp = processed_df[(1:top_ranking),3]
    processed_df[(1:top_ranking),3] = ifelse(temp==2,1,0)
    processed_df[,3][processed_df[,3]==2] = 0
    processed_df = arrange(processed_df,loc)
  }
  return(processed_df)
}

plot_df = function(processed_df,
                   lower_limit = min(processed_df$loc),
                   upper_limit = max(processed_df$loc),
                   peak_color = "coral",
                   peak_shape = 21,
                   peak_size = 1){
  processed_df = filter(processed_df,
                        loc>=lower_limit,
                        loc<=upper_limit)
  plt = ggplot(data = processed_df, aes(x = loc, y = score)) +
    geom_line() +
    geom_point(data = filter(processed_df,peak==1),
               color = peak_color,pch = peak_shape,
               size = peak_size) +
    xlab("Chromsome Location") + 
    ylab("Filtered Score") + 
    theme_gray() + 
    theme(text = element_text(size = 15))
  plt = ggplotly(plt)
  return(plt)
}

cms_list = read.csv("data.csv") %>%
  arrange(chr,loc) %>%
  select(chr,loc,cms_score)


