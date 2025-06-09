# Pacotes necessários ----
library(stringi)
library(Rcpp)

# Funções Auxiliares ----
## $N(w,a)$: Número de ocorrências de uma sequêcia seguida de um simbolo específico na amostra ----
N_w<- function(w, amostra, den = F) {
  if (nchar(w) == 0 || nchar(amostra) < nchar(w)) return(0)
  if(den == F) matches <- stri_locate_all_regex(paste0(amostra), paste0("(?=", w, ")"))[[1]]
  else matches <- stri_locate_all_regex(paste0(amostra), paste0("(?=", w, "\\d{1})"))[[1]]
  if (is.na(matches[1,1])) return(0) 
  else {
    return(nrow(matches))
  }
}

## $\hat{Q}(a|w)$: Probabilidade estimada ----
p_aw <- function(a, w, amostra){
  den <- N_w(w, amostra, den = T)
  if(den == 0){
    return(0)
  }
  else{
    prob <- N_w(paste0(w,a), amostra)/den
    return(prob)
  }
}

## Conjunto de todas as sequências finitas com tamanho até d que aparecem na amostra ----
calcula_Vn <- function(x, d) {
  s <- paste(x, collapse = "")
  Vn <- character()
  for (tam in 1:d) {
    substrings <- stri_sub(s, 1:(nchar(s) - tam + 1), length = tam)
    Vn <- union(Vn, substrings)
  }
  return(Vn)
}

### Conjunto de toda as sequencias finitas com tamanho até d ----
cataloga_d <- function(alfabeto, d) {
  if (d == 0) return("")
  do.call(expand.grid, rep(list(alfabeto), d)) |>
    apply(1, stri_paste, collapse = "")
}


## Sufixos de uma sequência ----
sufixos <- function(x, l_max = NULL, l_min = NULL, n = NULL) {
  if(is.null(n)) n <- nchar(x)
  if(is.null(l_max)) l_max <- n
  if(is.null(l_min)) l_min <- 1
  sufixos <- stri_sub(x, l_min:l_max, n)
  return(sufixos)
}

## Maior sufixo (pai) ----
pai <- function(x) {
  if(nchar(x) == 1) return('lambda')
  stri_sub(x, 2, nchar(x))
}

## Filhos de uma sequência ----
filhos <-  function(x, alfabeto) {
  if(x=='lambda') return(alfabeto)
  alfabeto %s+% x
}

# Irmãos ----
irmaos <- function(x, alfabeto) {
  filhos(pai(x), alfabeto)
}

# Calcula contagens aassociadas a determinados nós
contagens <- function(nos, alfabeto, amostra, folhas = NULL, lambda = F){
  n <- nchar(amostra)
  if (is.null(folhas)) {
    cont <- t(sapply(nos, \(w) sapply(alfabeto, \(a) N_w(paste0(w, a), amostra))))
    if(lambda) cont <- rbind('lambda' = sapply(alfabeto, function(a) N_w(a, amostra)), cont)
    return(cont)
  }
  
  tam_max <- max(nchar(nos))
  finais <- stri_sub(amostra, from = (n - tam_max + 1):(n - 1), to = n - 1)
  finais <- intersect(finais, nos)
  base <- union(folhas, finais)
  cont <- matrix(NA, nrow = length(nos), ncol = length(alfabeto), dimnames = list(nos, alfabeto))
  cont[base, ] <- t(sapply(base, \(w) sapply(alfabeto, \(a) N_w(paste0(w, a), amostra))))
  
  faltam <- setdiff(nos, base)
  for (w in rev(faltam)) {
    filhos <- paste0(w, alfabeto)
    cont[w, ] <- sapply(filhos, \(f) if (f %in% rownames(cont)) sum(cont[f, ]) else 0)
  }
  
  if(lambda) cont <- rbind('lambda' = sapply(alfabeto, function(a) N_w(a, amostra)), cont)
  return(cont)
}

## Simula de uma cadeia com memória de comprimento variável ----
cppFunction('
std::string simula(NumericMatrix transicoes, int N) {
  // Obtém os nomes das linhas (contextos)
  CharacterVector context_names = rownames(transicoes);
  int depth = 0;
  
  // Descobre profundidade máxima
  for (int i = 0; i < context_names.size(); i++) {
    int len = Rf_length(STRING_ELT(context_names, i));
    if (len > depth) depth = len;
  }

  // Descobre o alfabeto
  std::set<char> alpha_set;
  for (int i = 0; i < context_names.size(); i++) {
    std::string ctx = Rcpp::as<std::string>(context_names[i]);
    for (char c : ctx) alpha_set.insert(c);
  }

  std::vector<char> alphabet(alpha_set.begin(), alpha_set.end());
  int alph_size = alphabet.size();
  int burn_in = 1000 * alph_size;  // Tamanho do burn-in

  // Inicialização da saída (agora com espaço para o burn-in)
  std::vector<char> x(N + depth + burn_in);

  // Se profundidade 0, gera i.i.d.
  if (depth == 0) {
    NumericVector probs = transicoes(0, _);
    NumericVector cumsum_probs = Rcpp::cumsum(probs);
    
    // Gera burn-in
    for (int i = 0; i < burn_in; i++) {
      double u = R::runif(0, 1);
      for (int j = 0; j < alph_size; j++) {
        if (u <= cumsum_probs[j]) {
          x[i] = alphabet[j];
          break;
        }
      }
    }
    
    // Gera as observações principais
    for (int i = burn_in; i < N + burn_in; i++) {
      double u = R::runif(0, 1);
      for (int j = 0; j < alph_size; j++) {
        if (u <= cumsum_probs[j]) {
          x[i] = alphabet[j];
          break;
        }
      }
    }
  } else {
    // Início aleatório
    for (int i = 0; i < depth; i++) {
      x[i] = alphabet[rand() % alph_size];
    }

    // Fase de burn-in
    for (int i = depth; i < burn_in + depth; i++) {
      bool match_found = false;
      for (int d = depth; d >= 1; d--) {
        std::string ctx = "";
        for (int j = d; j >= 1; j--) {
          ctx += x[i - j];
        }

        // Procura o contexto na matriz
        for (int k = 0; k < context_names.size(); k++) {
          if (Rcpp::as<std::string>(context_names[k]) == ctx) {
            NumericVector probs = transicoes(k, _);
            NumericVector cumsum_probs = Rcpp::cumsum(probs);
            double u = R::runif(0, 1);
            for (int j = 0; j < alph_size; j++) {
              if (u <= cumsum_probs[j]) {
                x[i] = alphabet[j];
                break;
              }
            }
            match_found = true;
            break;
          }
        }
        if (match_found) break;
      }
      if (!match_found) {
        stop("Contexto não encontrado na posição " + std::to_string(i));
      }
    }

    // Gera as observações principais após o burn-in
    for (int i = burn_in + depth; i < N + burn_in + depth; i++) {
      bool match_found = false;
      for (int d = depth; d >= 1; d--) {
        std::string ctx = "";
        for (int j = d; j >= 1; j--) {
          ctx += x[i - j];
        }

        // Procura o contexto na matriz
        for (int k = 0; k < context_names.size(); k++) {
          if (Rcpp::as<std::string>(context_names[k]) == ctx) {
            NumericVector probs = transicoes(k, _);
            NumericVector cumsum_probs = Rcpp::cumsum(probs);
            double u = R::runif(0, 1);
            for (int j = 0; j < alph_size; j++) {
              if (u <= cumsum_probs[j]) {
                x[i] = alphabet[j];
                break;
              }
            }
            match_found = true;
            break;
          }
        }
        if (match_found) break;
      }
      if (!match_found) {
        stop("Contexto não encontrado na posição " + std::to_string(i));
      }
    }
  }

  // Concatena o output (ignorando o burn-in e os valores iniciais de depth)
  std::string output = "";
  int start_index = burn_in + depth;
  for (int i = start_index; i < N + start_index; i++) {
    output += x[i];
  }

  return output;
}
')







