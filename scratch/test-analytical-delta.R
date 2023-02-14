library(tidyverse)
source("dist_structural_util.R")
library(microbenchmark)
library(truncnorm)
library(testthat)
library(cmdstanr)
library(rstan)
library(posterior)
library(tidybayes)
library(nleqslv)




expose_cmdstanr_functions <- function(model_path, include_paths = NULL,
                                    expose_to_global_env = FALSE) {
  required_pkgs <- c("Rcpp", "RcppEigen", "cmdstanr")
  found_pkgs <- required_pkgs %in% rownames(installed.packages())
  if (!all(found_pkgs)) {
    stop(
      "The following required packages are missing: ",
      paste0(required_packages[!found_pkgs], collapse = ", "),
      "."
    )
  }
  if (cmdstanr::cmdstan_version() < "2.26.0") {
    stop("Please install CmdStan version 2.26 or newer.", call. = FALSE)
  }
  get_cmdstan_flags <- function(flag_name) {
    cmdstan_path <- cmdstanr::cmdstan_path()
    flags <- processx::run(
      "make", 
      args = c(paste0("print-", flag_name)),
      wd = cmdstan_path
    )$stdout
    flags <- gsub(
      pattern = paste0(flag_name, " ="),
      replacement = "", x = flags, fixed = TRUE
    )
    flags <- gsub(
      pattern = " stan/", replacement = paste0(" ", cmdstan_path, "/stan/"),
      x = flags, fixed = TRUE
    )
    flags <- gsub(
      pattern = "-I lib/", replacement = paste0("-I ", cmdstan_path, "/lib/"),
      x = flags, fixed = TRUE
    )
    flags <- gsub(
      pattern = "-I src", replacement = paste0("-I ", cmdstan_path, "/src"),
      x = flags, fixed = TRUE
    )
    gsub("\n", "", flags)
  }
  temp_stan_file <- tempfile(pattern = "model-", fileext = ".stan")
  temp_cpp_file <- paste0(tools::file_path_sans_ext(temp_stan_file), ".cpp")
  file.copy(model_path, temp_stan_file, overwrite = TRUE)
  if (isTRUE(.Platform$OS.type == "windows")) {
    stanc3 <- "./bin/stanc.exe"
  } else {
    stanc3 <- "./bin/stanc"
  }
  processx::run(
    stanc3,
    args = c(
      temp_stan_file,
      "--standalone-functions",
      paste0("--include-paths=", include_paths),
      paste0("--o=",temp_cpp_file)
    ),
    wd = cmdstanr::cmdstan_path()
  )
  code <- paste(readLines(temp_cpp_file), collapse = "\n")
  code <- paste(
    "// [[Rcpp::depends(RcppEigen)]]",
    "#include <stan/math/prim/fun/Eigen.hpp>",
    "#include <RcppCommon.h>
    #include <boost/random/additive_combine.hpp>
    #include <iostream>

    namespace Rcpp {
      SEXP wrap(boost::ecuyer1988 RNG);
      SEXP wrap(boost::ecuyer1988& RNG);
      SEXP wrap(std::ostream stream);
      template <> boost::ecuyer1988 as(SEXP ptr_RNG);
      template <> boost::ecuyer1988& as(SEXP ptr_RNG);
      template <> std::ostream* as(SEXP ptr_stream);
      namespace traits {
        template <> class Exporter<boost::ecuyer1988&>;
        template <> struct input_parameter<boost::ecuyer1988&>;
      }
    }

    #include <Rcpp.h>

    namespace Rcpp {
      SEXP wrap(boost::ecuyer1988 RNG){
        boost::ecuyer1988* ptr_RNG = &RNG;
        Rcpp::XPtr<boost::ecuyer1988> Xptr_RNG(ptr_RNG);
        return Xptr_RNG;
      }

      SEXP wrap(boost::ecuyer1988& RNG){
        boost::ecuyer1988* ptr_RNG = &RNG;
        Rcpp::XPtr<boost::ecuyer1988> Xptr_RNG(ptr_RNG);
        return Xptr_RNG;
      }

      SEXP wrap(std::ostream stream) {
        std::ostream* ptr_stream = &stream;
        Rcpp::XPtr<std::ostream> Xptr_stream(ptr_stream);
        return Xptr_stream;
      }

      template <> boost::ecuyer1988 as(SEXP ptr_RNG) {
        Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
        boost::ecuyer1988& RNG = *ptr;
        return RNG;
      }

      template <> boost::ecuyer1988& as(SEXP ptr_RNG) {
        Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
        boost::ecuyer1988& RNG = *ptr;
        return RNG;
      }

      template <> std::ostream* as(SEXP ptr_stream) {
        Rcpp::XPtr<std::ostream> ptr(ptr_stream);
        return ptr;
      }

      namespace traits {
        template <> class Exporter<boost::ecuyer1988&> {
        public:
          Exporter( SEXP x ) : t(Rcpp::as<boost::ecuyer1988&>(x)) {}
          inline boost::ecuyer1988& get() { return t ; }
        private:
          boost::ecuyer1988& t ;
        } ;

        template <>
        struct input_parameter<boost::ecuyer1988&> {
          typedef
          typename Rcpp::ConstReferenceInputParameter<boost::ecuyer1988&> type ;
          //typedef typename boost::ecuyer1988& type ;
        };
      }
    }

    RcppExport SEXP get_stream_() {
      std::ostream* pstream(&Rcpp::Rcout);
      Rcpp::XPtr<std::ostream> ptr(pstream, false);
      return ptr;
    }

    RcppExport SEXP get_rng_(SEXP seed) {
      int seed_ = Rcpp::as<int>(seed);
      boost::ecuyer1988* rng = new boost::ecuyer1988(seed_);
      Rcpp::XPtr<boost::ecuyer1988> ptr(rng, true);
      return ptr;
    }
    ",
    "#include <RcppEigen.h>",
    code,
    sep = "\n"
  )
  code <- gsub("// [[stan::function]]",
              "// [[Rcpp::export]]", code, fixed = TRUE)
  code <- gsub(
    "stan::math::accumulator<double>& lp_accum__, std::ostream* pstream__ = nullptr){",
    "std::ostream* pstream__ = nullptr){\nstan::math::accumulator<double> lp_accum__;",
    code,
    fixed = TRUE
  )
  code <- gsub("__ = nullptr", "__ = 0", code, fixed = TRUE)

  get_stream <- function() {
    return(.Call('get_stream_'))
  }
  get_rng <- function(seed=0L) {
    if (!identical(seed, 0L)) {
      if (length(seed) != 1)
        stop("Seed must be a length-1 integer vector.")
    }
    return(.Call('get_rng_', seed))
  }
  if (expose_to_global_env) {
    env = globalenv()
  } else {
    env = new.env()
  }
  compiled <- withr::with_makevars(
    c(
      USE_CXX14 = 1,
      PKG_CPPFLAGS = "",
      PKG_CXXFLAGS = get_cmdstan_flags("CXXFLAGS"),
      PKG_LIBS = paste0(
        get_cmdstan_flags("LDLIBS"),
        get_cmdstan_flags("LIBSUNDIALS"),
        get_cmdstan_flags("TBB_TARGETS"),
        get_cmdstan_flags("LDFLAGS_TBB")
      )
    ),
    Rcpp::sourceCpp(code = code, env = env)
  )
  for (x in compiled$functions) {
    FUN <- get(x, envir = env)
    args <- formals(FUN)
    args$pstream__ <- get_stream()
    if ("lp__" %in% names(args)) args$lp__ <- 0
    if ("base_rng__" %in% names(args)) args$base_rng__ <- get_rng()
    formals(FUN) <- args
    assign(x, FUN, envir = env)
  }
  assign("stan_rng__", get_rng, envir = env)
  if (expose_to_global_env) {
    invisible(NULL)
  } else {
    return(env)
  }
}



exposed_funcs = expose_cmdstanr_functions(model_path = "stan_models/ed-stan-funcs.stan") 



r_data = read_rds(
  "temp-data/r-predicted-data.rds"
)

source("optim/optim-functions.R")

r_pred_df = r_data$pred_df %>%
  mutate(
    u_sd = sqrt(total_error_sd^2 - 1)
  )

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)
# stick to monitored sms.treatment group
# remove sms.treatment.2
monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data %>% 
mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group,
    sms_treatment = factor(sms.treatment.2))


# fit_file = "data/stan_analysis_data/dist_fit71_STRUCTURAL_LINEAR_U_SHOCKS-1.csv"

# fit = as_cmdstan_fit(fit_file)


# To find fixed point need:
# benefit_cost
# mu_rep
# total_error_sd
# u_sd



# bc_draws = fit$draws(
#   variables = c(
#     "structural_cluster_benefit_cost", 
#     "cluster_dist_cost", 
#     "obs_cluster_mu_rep",
#     "total_error_sd[1]", 
#     "u_sd[1]")
# )


# bc_draws %>%
#   saveRDS("temp-data/temp-bc-draws.rds")

bc_draws = read_rds("temp-data/temp-bc-draws.rds")   
  

# rm(fit)
# gc()
# stop()



bc_draw_df = bc_draws %>%
  as_draws_df()

bc_rvar_df = bc_draw_df %>%
  spread_rvars(
    structural_cluster_benefit_cost[cluster_id], 
    cluster_dist_cost[cluster_id],
    obs_cluster_mu_rep[cluster_id],
    total_error_sd, 
    u_sd
    ) %>%
    left_join(
      analysis_data %>%
          select(
            cluster_id,
            assigned_dist_group, 
            assigned_treatment = assigned.treatment, 
            dist = cluster.dist.to.pot) %>%
            unique(), 
          by = "cluster_id"
    )  


draw_bc_rvar_df = bc_rvar_df %>%
  unnest_rvars()

draw_bc_rvar_df
library(furrr)
plan(multicore, workers = 8)
draw_bc_rvar_df = draw_bc_rvar_df %>%
  mutate(
      R_v_star = future_pmap(
        list(
          dist, 
          structural_cluster_benefit_cost,
          obs_cluster_mu_rep, 
          total_error_sd, 
          u_sd
        ),
        ~find_v_star(
          distance = ..1, 
          b = ..2, 
          mu_rep = ..3, 
          total_error_sd = ..4, 
          u_sd = ..5, 
          bounds = c(-Inf, Inf)
      )
  )
)

draw_bc_rvar_df = draw_bc_rvar_df %>%
  mutate(
      stan_v_star = future_pmap_dbl(
        list(
          structural_cluster_benefit_cost,
          obs_cluster_mu_rep, 
          total_error_sd, 
          u_sd
        ),
        ~exposed_funcs$find_fixedpoint_solution(
          b = ..1, 
          mu_rep = ..2, 
          total_error_sd = ..3, 
          u_sd = ..4, 

          use_u_in_delta = 1,
          alg_sol_f_tol = 0.001,
          alg_sol_max_steps = 1e9L,
          alg_sol_rel_tol = 0.0000001

      ),
      .options = furrr_options(seed = TRUE), 
      .progress = TRUE
    )
  )
draw_bc_rvar_df = draw_bc_rvar_df %>%
  unnest_wider(R_v_star)


draw_bc_rvar_df






draw_bc_rvar_df = draw_bc_rvar_df %>%
  mutate(
    pred_takeup = 1 - pnorm(v_star/total_error_sd)
  )



draw_bc_rvar_df %>%
  group_by(
    cluster_id, 
    dist, 
    assigned_treatment
  ) %>%
  summarise(
    v_star = mean(v_star), 
    takeup = mean(pred_takeup), 
    stan_takeup = mean(1 - pnorm(stan_v_star/total_error_sd)),
    stan_v_star = mean(stan_v_star)
  ) %>%
  select(-contains("takeup")) %>%
  gather(
    variable, 
    value, 
    v_star, stan_v_star
    ) %>%
  ggplot(aes(
    x = dist, 
    y = value, 
    colour = variable
  )) +
  geom_point() +
  geom_line() +
  facet_wrap(~assigned_treatment)  +
  labs(
    title = "Using Parameter Posterior Draws", 
    subtitle = "Fit using exposed Stan funcs and R"
  )


ggsave(
  "temp-plots/post-draws-comp.png",
  width = 8, 
  height = 6, 
  dpi = 500
)


clust_bs_df = bc_draw_df %>%
  as_tibble() %>%
  gather(variable, value) %>%
  filter(!str_detect(variable, "\\.")) %>%
  group_by(variable) %>% 
  mutate(draw = 1:n()) %>%
  mutate(
    cluster_id = str_extract(variable, "\\d+") %>% as.integer
  )  %>%
  left_join(
    analysis_data %>%
        select(cluster_id,assigned_dist_group, assigned.treatment, cluster.dist.to.pot), 
        by = "cluster_id"
  )  %>%
  mutate(
    variable = str_remove(variable, "\\[.*$")
  )







b_stan_df = clust_bs_df %>%
  group_by(cluster_id, variable) %>%
  summarise(
    mean_est = mean(value), 
    treatment = unique(assigned.treatment), 
    dist =  unique(cluster.dist.to.pot)
  )


  
 b_stan_df %>% 
  ggplot(aes(
    x = dist, 
    y = mean_est, 
    colour =treatment 
  )) +
  geom_point() +
  facet_wrap(~variable)

b_stan_df




b_comp_df = bind_rows(
  r_pred_df  %>%
    select(dist, b, treatment, draw) %>%
    group_by(dist, treatment) %>%
    summarise(
      mean_est = mean(b), 
      median_est = median(b)
    ) %>%
    mutate(type = "ed"), 
    b_stan_df %>% 
      filter(variable == "structural_cluster_benefit_cost") %>%
      mutate(type = "stan")
)





r_data$param_grid %>%
  select(k) %>%
  unique()

r_data$param_grid %>%
  filter(.variable == "beta") %>%
  group_by(k) %>%
  summarise(
    median_draw = median(.value)
  )

  r_data$param_grid %>%
    filter(.variable == "beta") %>%
      filter(k == 2) %>%
      select(draw = .draw, beta_control = .value)

b_comp_df %>%
  ggplot(aes(
    x = dist, 
    y = mean_est, 
    colour = type
  )) +
  facet_wrap(~treatment) +
  geom_point() +
  geom_line() +
  labs(
    title = "Net Private Benefit Stan Posterior Draws vs R function of Post draws"
  )

ggsave(
  "temp-plots/temp.png", 
  width = 8, 
  height = 6, 
  dpi = 500
)

b_comp_df %>%
  ggplot(aes(
    x = dist, 
    y = mean_est, 
    colour = treatment
  )) +
  facet_wrap(~type) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom")





r_pred_df %>%
  group_by(
    dist, treatment
  ) %>%
    summarise(
      b = mean(b)
    )

































draws_used = r_data$pred_df$draw %>% unique()


r_pred_df = r_data$pred_df %>%
  mutate(
    u_sd = sqrt(total_error_sd^2 - 1)
  )

r_pred_df = r_pred_df %>%
  mutate(
    expected_delta = pmap_dbl(
      list(
        v_star, 
        total_error_sd, 
        u_sd
      ), 
    ~exposed_funcs$expected_delta(
      w = ..1, 
      total_error_sd = ..2, 
      u_sd = ..3, 
      x_r = 1, 
      x_i = 1
    )
    ) 
  )


r_pred_df

r_pred_df %>%
  ggplot(aes(
    x = expected_delta, 
    y = delta_v_star
  )) +
  geom_point() +
  geom_abline(
    linetype = "longdash"
  )


poss_stan_fp = possibly(
  exposed_funcs$find_fixedpoint_solution,
  otherwise = NA 
)
stan_fp = exposed_funcs$find_fixedpoint_solution

r_pred_df = r_pred_df %>%
  mutate(
      stan_v_star = pmap_dbl(
        list(
          b, 
          mu_rep,
          total_error_sd, 
          u_sd
        ), 
        ~stan_fp(
          benefit_cost = ..1, 
          mu_rep = ..2, 
          total_error_sd = ..3, 
          u_sd = ..4, 

          use_u_in_delta = 1,
          alg_sol_f_tol = 0.001,
          alg_sol_max_steps = 1e9L,
          alg_sol_rel_tol = 0.0000001
        )
      )
  )


r_pred_df %>%
  filter(is.na(stan_v_star)) %>%
  select(v_star, delta_v_star)


r_pred_df = r_pred_df %>%
  mutate(
    stan_pred_takeup = 1 - pnorm(stan_v_star / total_error_sd)
  )


r_pred_df %>%
  select(dist, treatment, stan_v_star, v_star, draw) %>%
  gather(variable, value, stan_v_star, v_star)   %>%
  group_by(
    variable, treatment, dist
  ) %>%
  summarise(
    mean_est = mean(value), 
    median_est = median(value)
  ) %>%
  ggplot(aes(x = dist, y = mean_est, colour = variable)) +
  geom_point() +
  facet_grid(variable ~ treatment)

r_pred_df %>%
  select(dist, treatment, stan_v_star, v_star, draw) %>%
  gather(variable, value, stan_v_star, v_star)  %>%
  # filter(draw %in% 1:200) %>%
  ggplot(aes(
    x = dist, 
    y = value, 
    colour = variable
  )) +
  facet_grid(variable~treatment) +
  geom_point()




r_pred_df %>%
  select(dist, treatment, stan_pred_takeup, pred_takeup, draw) %>%
  gather(variable, value, stan_pred_takeup, pred_takeup)  %>%
  filter(draw == sample(r_data$pred_df$draw, 1)) %>%
  ggplot(aes(
    x = dist, 
    y = value, 
    colour = variable
  )) +
  facet_wrap(~treatment) +
  geom_point()


r_pred_df %>%
  ggplot(aes(
    x = pred_takeup, 
    y = stan_pred_takeup, 
    colour = treatment
  )) +
  geom_point() +
  geom_abline(linetype = "longdash")

r_pred_df %>%
  ggplot(aes(
    x = v_star, 
    y = stan_v_star
  )) +
  geom_point() +
  geom_abline(linetype = "longdash")

















stan_owen_t_code <-
  '
  functions {
    vector stan_owen_t(vector x, vector y) {
      return owens_t(x, y);
   }
  }
'
expose_stan_functions(stanc(model_code = stan_owen_t_code))


test_that("Delta Equal at limits", {
  ad = analytical_delta(3, 0.1)
  ad_below = analytical_delta_bounded(3, 0.1, bounds = c(-3, 3))
  ad_b = analytical_delta_bounded(3, 0.1, bounds = c(-10000, 10000))

  ad_b_inf = analytical_delta_bounded(3, 0.1, bounds = c(-Inf, Inf))


  expect_equal(
    ad, ad_b
  )

  expect_equal(
    ad,
    ad_b_inf
  )

  expect_lte(ad_below, 3)

})



ub = 3
lb = -3
bandwidth = 3
df = expand.grid(
  w = seq(from = lb - bandwidth, to = ub + bandwidth, length.out = 100 ),
  u_sd = seq(from = 0.1, to = 2,  by = 0.05)
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    delta = calculate_delta(w, sqrt(1 + u_sd^2), u_sd),
    analytical_delta = analytical_delta(w, u_sd), 
    delta_deriv = calculate_delta_deriv(w, sqrt(1 + u_sd^2), u_sd),
    analytical_delta_deriv = analytical_delta_deriv(w, u_sd, delta_w = analytical_delta),
    delta_bounded = calculate_delta_bounded(w, sqrt(1 + u_sd^2), u_sd, bounds = c(lb, ub)),
    analytical_delta_bounded = analytical_delta_bounded(w, u_sd, bounds = c(lb, ub)),
    analytical_delta_deriv_bounded = analytical_delta_deriv_bounded(w, u_sd, bounds = c(lb, ub), delta_w = analytical_delta_bounded),
    analytical_conv_Fw = analytical_conv_Fw(w, u_sd, bounds = c(lb, ub))
  )




comp_df = df %>%
  pivot_longer(
    delta:analytical_conv_Fw
  ) %>%
  mutate(
    type = if_else(str_detect(name, "analytical_"), "analytical", "numerical"), 
    name = str_remove(name, "analytical_")
  ) %>%
  spread(type, value)

make_plot = TRUE
if (make_plot == TRUE) {
  p_comp_plot = comp_df %>%
    filter(!(name %in% c("delta_deriv_bounded", "conv_Fw"))) %>%
    ggplot(aes(
      x = analytical, 
      y = numerical, 
      colour = name
    )) +
    geom_point() +
    facet_wrap(~name, scales = "free", ncol = 1) +
    theme_bw() +
    guides(colour = "none") +
    geom_abline(linetype = "longdash") +
    labs(
      title = "Numerical vs Analytical Delta Calculations", 
      subtitle = str_glue("Delta bounded by {lb}, {ub}")
    )
  p_comp_plot
  ggsave(
    plot = p_comp_plot,
    "temp-plots/numerical-comp-plot.png",
    width = 8,
    height = 6,
    dpi = 500
  )

  comp_df %>%
    gather(variable, value, analytical, numerical ) %>%
    filter(variable == "analytical") %>%
    filter(name == "delta_deriv_bounded") %>%
    ggplot(aes(
      x = w, 
      y = pmin(value, 3), 
      colour = variable, 
      group = u_sd
    )) +
    facet_wrap(~variable) + 
    geom_line()  +
    geom_vline(xintercept = c(lb, ub), linetype = "longdash") +
    guides(colour = "none") +
    theme_bw() +
    labs(
      title = "Derivative of Delta, Bounded"
    )
  ggsave("temp-plots/delta-deriv-bounded.png", width = 8,  height = 6, dpi = 500)

  comp_df %>%
    gather(variable, value, analytical, numerical ) %>%
    filter(name == "delta_bounded") %>%
    ggplot(aes(
      x = w, 
      y = pmin(value, 3), 
      colour = variable, 
      group = u_sd
    )) +
    facet_wrap(~variable) + 
    geom_line()  +
    theme_bw() +
    guides(colour = "none") +
    geom_vline(xintercept = c(lb, ub), linetype = "longdash") +
    labs(title = "Comparing Delta, Bounded")

ggsave("temp-plots/delta-bounded-numerical-comp.png",
  width = 8,
  height = 6,
  dpi = 500)

}


summ_comp_df = comp_df %>%
  filter(name != "delta_deriv_bounded") %>%
  mutate(
    diff = analytical - numerical
  ) %>%
  group_by(name) %>%
  summarise(
    median_diff = median(diff), 
    mean_diff = mean(diff), 
    median_abs_diff = median(abs(diff)), 
    mean_abs_diff = mean(abs(diff))
  )

test_that("Analytical and numerical align", {
  mad = summ_comp_df$median_abs_diff
  mean_abs_diff = summ_comp_df$mean_abs_diff
  map(
    mad,
    expect_lte, 1e-4
  )
  # Deriv a bit different
  map(
    mean_abs_diff,
    expect_lte, 1e-1
  )
})


benchmark = FALSE
if (benchmark == TRUE) {
  benchmark_results = microbenchmark(
    "analytical_delta" = analytical_delta(seq(from = -3, to = 3, by = 0.1), 0.2), 
    "numerical_delta" = calculate_delta(seq(from = -3, to = 3, by = 0.1), sqrt(1 + 0.2^2), 0.2), 
    "analytical_delta_deriv" = analytical_delta_deriv(seq(from = -3, to = 3, by = 0.1), 0.2), 
    "numerical_delta_deriv" = calculate_delta(seq(from = -3, to = 3, by = 0.1), sqrt(1 + 0.2^2), 0.2),
    "analytical_delta_bounded" = 
        analytical_delta_bounded(seq(from = -3, to = 3, by = 0.1),
          0.2, c(-3, 3)),
    "numerical_delta_bounded" = calculate_delta_bounded(seq(from = -3, to = 3, by = 0.1), sqrt(1 + 0.2^2), 0.2, c(-3, 3)),
    "analytical_delta_deriv_bounded" = 
        analytical_delta_deriv_bounded(seq(from = -3, to = 3, by = 0.1),
          0.2, c(-3, 3))
  )

  benchmark_results
}



