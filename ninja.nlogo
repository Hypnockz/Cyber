;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Versión de EIA 2.0(AB), simplificado basado en topología de AC.
;;Proyecto financiado por Universidad Santo Tomás a través de fondo de Investigación y Creación
;;Trabajo realizado en el contexto de Tesis de Ph.D. Bajo la Supervisión de los Doctores Lozano y Blum de la EHU, España.
;;Derechos reservados a Pedro Pinacho Davidson 2013.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Observaciones
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Esta versión contempla, las siguientes variantes:
;;-Heurística para castigar vecinos similares.
;;-Todos se alimentan del mismo alimento.
;;-Reducción de dimensión a 1D
;;-Posibilidad de uso de envejecimiento y muerte
;;-Versión para usar procesamiento batch de archivos
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Cambio en esta versión
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;-Los agentes son cambiados por celdas con estado determinado por energía y vector genético asociado.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;------------------------------------------------------------------------------------------------------------------------
;;Descripción de variables globales
;;------------------------------------------------------------------------------------------------------------------------
;;parámetros del modelo previamente ajustados mediante metaaprendizaje con AG tutor.

;;energia-inicial: parametro del model que establece la energía inicial de los agentes.
;;maximo-aporte-energia-sedimento: parámetro del modelo que establece la máxima energía entregada por la partícula alpha en caso de mayor afinidad.
;;toxicidad-incremental: parámetro del modelo que establece la pérdida  lineal del aporte de energía por falta de afinidad genética del agente hacia la partícula.
;;consumo-energia-basal: descuento de energía a todo agente por iteración
;;razon-reproduccion: diferencia que debe haber de energía para transformar el estado de una celda por el del vecino
;;penalizacion-similitud: penalización de energía a celdas similares
;;particula: no utilizado en esta versión
;;energia-maxima:  energía maxima evitando insensibilizar modelo
;;particulas: contenedor de partículas alpha
;;longitud-automata: longitud de autómata  lineal que representa a los agentes

;;-------------------------------------------------------------------------------------------------------------------------
;;variables del clasificador
;;energia-promedio-algas: variable para propósitos gráficos y para uso en proceso de clasificación
;;variable-ventana-deslizante-energia: arreglo que almacena energías promedio de la población de agentes implementando ventana deslizante
;;variable-ventana-deslizante-poblacion: arreglo que almacena tamaño de la población de agentes implementando ventana deslizante
;;media-muestral-deslizante-energia:  media de ventana deslizante de energía promedio de población
;;media-muestral-deslizante-poblacion:  media de ventana deslizante  de tamaño de población
;;desviacion-muestral-deslizante-energia:   desviación de ventana deslizante de energía promedio de población
;;desviacion-muestral-deslizante-poblacion: desviación de ventana deslizante  de tamaño de población
;;bandas-xxx-xxx: bandas límite que bordean las tendencias de la ventana deslizante, usadas para detectar anomalías a la tendencia.
;;ataque-detectado: indica traspaso de umbral de tolerancia, por tiempo acumulado traspasando límites de bandas (salida-banda) (ver gráfico alertas tiempo)
;;salida-banda: contador de iteraciones consecutivas donde las bandas de las ventanas deslizantes son traspasadas por la media de energía o tamaño de población.
;;umbral-alerta: parámetro de clasificador que indica el valor que debe obtener salida-banda para volver cierto ataque-detectado
;;VP, FP, FN, VN: Verdadero Positivo, Falso Positivo, Falso Negativo y Verdadero Negativo respectivamente.
;;sensibilidad, especificidad: indicadores calculados a partir de tabla de contingencia.


;;-------------------------------------------------------------------------------------------------------------------------
;;otras variables
;;ataque: variable usada para propósitos gráficos.
;;archivos-procesados: cantidad de archivos de tráfico procesados.
;;caracteristicas-archivo: arreglo contenedor de características cargadas desde archivo de datos
;;indice-caracteristicas-utilizadas: arreglo con indicadores de características binarias del tráfico utilizadas por el clasificador
;;fuente-datos, fuente-datos-bin: variables para salida a usuario y señalización de tipos de datos analizados
;;lapso-muestreo: indica granularidad utilizada para el cálculo de tabla de contingencia del clasificador (unitario en esta versión)
;;muestras-revisadas: variable utilizada para brindar información al usuario
;;contador-iteraciones: no necesita explicación
;;entrenamiento: variable utilizada para dar información al usuario
;;clasificacion: variable accesoria utilizada para la síntesis de la tabla de contingencia.
;;flag-par nombre-archivo conteo-archivo: variables accesorias usadas para control de archivos de datos.
;;numeros-primos:  lista de números primos utlizados con propósitos de cálculo de variabilidad-genética
;;profundidad-histograma: variables propias del histograma de variabilidad genética
;;indice-diversidad: variables propias del histograma de variabilidad genética
;;numero-genotipos: variables propias del histograma de variabilidad genética
;;flag-diversidad: variables propias del histograma de variabilidad genética
;;media-normal: no utilizado en esta versión
;;desviacion-normal: no utilizado en esta versión
;;energia-acumulado:  variable con energía acumulada, usada para propósitos gráficos
;;celdas-amarillas : cantidad de celdas en estado amarillo, usada para propósitos gráficos
;;celdas-rojas: cantidad de celdas en estado crítico de energía, usada para propósitos gráficos
;;celdas-verdes: cantidad de celdas en estado de alta energía, usada para propósitos gráficos
;;celdas-totales: cantidad total de celdas (vivas)
;;limite-energia: límite de banda de energía para evitar insensibilización del model (experimental)
;;limite-poblacion: límite de banda de tamaño población para evitar insensibilización del model (experimental)
;;factor-aumento-alerta: factor de incremento de alerta por salida de bandas
;;factor-disminucion-alerta: factor de disminución de alerta por retorno a la banda
;;fin-simulacion: variable de control de simulación.

Globals [energia-inicial maximo-aporte-energia-sedimento toxicidad-incremental consumo-energia-basal
  energia-promedio-algas estado-sedimentos rango-sedimentos variable-histograma numeros-primos
  profundidad-histograma indice-diversidad numero-genotipos flag-diversidad media-normal
  desviacion-normal variable-ventana-deslizante-energia variable-ventana-deslizante-poblacion
  energia-nuevo-agente poblaciones indicador-variabilidad media-muestral-deslizante-energia
  media-muestral-deslizante-poblacion desviacion-muestral-deslizante-energia
  desviacion-muestral-deslizante-poblacion banda-superior-energia banda-inferior-energia
  alerta-anomalia-energia banda-superior-poblacion banda-inferior-poblacion
  alerta-anomalia-poblacion ataque-detectado salida-banda umbral-alerta
  energia-acumulado

;; Nuevas Variables

razon-reproduccion penalizacion-similitud  particula energia-maxima particulas
longitud-automata
;; Gráficos

celdas-amarillas celdas-rojas  celdas-verdes celdas-totales

;; Manejo de archivos con datos de ataques/normalidad
ataque VP FP FN VN clasificacion sensibilidad especificidad flag-par nombre-archivo muestras-revisadas entrenamiento contador-alertas fuente-datos fuente-datos-bin archivos-procesados
caracteristicas-archivo indice-caracteristicas-utilizadas
nuevo tipo lineas-leidas

PAR-banda-superior-energia PAR-banda-inferior-energia PAR-banda-superior-poblacion PAR-banda-inferior-poblacion
PAR-cantidad-celdas ventana-entrenamiento

;;EXPERIMENTALES
;;limites de tolerancia luego de adaptacion
limite-energia
limite-poblacion
factor-aumento-alerta
factor-disminucion-alerta
fin-simulacion


  ;;Agregadas paara versión con LevelSpace
  rounds
  activeDelay1
  activeDelay2

  Test1
  pob-prev
  pob-delta
  pheromones
  cryogenic-capsule
  flag_capsule
  pher-prom
  Tipo_ataque_value
  List_ss
  index_attack
  energy_prom
  energy_umbral
  energy_actual
  timer_attack
  timer_normality
  under_attack
 ]

;; celdas tienen energía asociada y vector genético
patches-own [energia genes]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Rutinas de inicialización
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to setup
  file-close-all
  clear-all
  set-patch-size 2
  ;set max-pxcor
  ask patches [set pcolor white]
  iniciar-parametros
  iniciar-celdas
  reset-ticks
end





to iniciar-parametros
  set archivos-procesados 0
  set ataque 0
  set energia-inicial 200                ;; energia "primeros agentes"
  set media-normal 12                    ;; media inicial función generadora
  set desviacion-normal 0.2              ;; desviación de función generadora
  set rango-sedimentos 25                ;; dispersión de la función generadora
  set estado-sedimentos 1                ;; Posición inicial de ventana de sedimentos
  set maximo-aporte-energia-sedimento 200 ;; Maximo entregado con afinidad total
  set consumo-energia-basal 1            ;; consumo de energía metabólico por iteración
  set toxicidad-incremental  8.8 ;5.5          ;; Unidad de incremento de toxicidad nutricional
  set razon-reproduccion 0.243061             ;; razón entre celda local  y vecindad para cambiar por influencia de la otra celda
  set energia-nuevo-agente 40            ;; energía inicial de celda recién cambiada por influencia de otra celda
  set energia-maxima 800                 ;; energía maxima evitando insensibilizar modelo
  set celdas-verdes 0                    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  set celdas-amarillas 0                 ;; Inicialización de gráficos
  set celdas-rojas 0                     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  set longitud-automata 60
  set penalizacion-similitud 100
  set Test1 69
  set Tipo_ataque_value 0
  set pob-prev  0
  set pob-delta  0
  set pheromones 0
  set cryogenic-capsule []
  set flag_capsule 0
  set pher-prom []
  set List_ss [3906 11276 13516 24053 42895 47190 50392 53060 55934 58180 62381 66405 69666 72447 75490 77508 80265 83190 86985 89178]
  set index_attack 0
  set energy_prom []
  set energy_umbral 0 ;; energy prom anterior
  set energy_actual 0
  set timer_attack 0
  set timer_normality 0
  set under_attack 0
  ;;Parámetros histograma


  set numeros-primos [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131
    137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 ]
  set profundidad-histograma 10
  set variable-histograma []


  ;;Preparación de ventanas deslizantes

  set variable-ventana-deslizante-energia [ 480 ]
  repeat 500 [set variable-ventana-deslizante-energia lput 480 variable-ventana-deslizante-energia]


  set variable-ventana-deslizante-poblacion [ 100 ]
  repeat 500 [set variable-ventana-deslizante-poblacion lput 100 variable-ventana-deslizante-poblacion]


  ;; Parámetros para implementación de clasificador

  set PAR-banda-superior-energia 100
  set PAR-banda-inferior-energia 94
  set PAR-banda-superior-poblacion 100
  set PAR-banda-inferior-poblacion 9
  set PAR-cantidad-celdas 100

  set energia-acumulado 0
  set umbral-alerta 53
  set ventana-entrenamiento 200

  set limite-energia 0
  set limite-poblacion 0
  set factor-aumento-alerta 7
  set factor-disminucion-alerta 8

  ;;Parámetros para uso de archivos de datos
  set entrenamiento "Entrenando"
  set contador-alertas 0
  set clasificacion -1
  set fuente-datos "no detectados"
  set fuente-datos-bin -1
  ;; Par�metros  para cargar datos desde archivo
  file-close
  file-open "datosNormalesEIA3000.dat"
  set fuente-datos "NORMALES ESTACIÓN"
  set fuente-datos-bin 0
  set VN 0
  set VP 0
  set FN 0
  set FP 0
  set sensibilidad 0
  set especificidad 0
  set muestras-revisadas 0

  set particulas []
  set lineas-leidas 0
  set fin-simulacion false

  set rounds 0
  set activeDelay1 0
  set activeDelay2 0

end





  to iniciar-celdas

    ask patches with [pxcor = 0] [  set genes  shuffle [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
         22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53]
    set energia energia-inicial
    set pcolor green
     ;;set pcolor energia
    ]
  end

to delta-poblacion

  let  pob-actual count patches with [pcolor != white and pxcor = 0]
  ;show pob-actual
  ;show pob-prev
  ;show pob-delta
  ;show "--"
  set pob-delta (pob-actual - pob-prev)

end




  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Rutina de Ejecución Principal
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go
  ifelse (not fin-simulacion)[
    if (ticks =  ventana-entrenamiento) [establecer-limites-bandas]
    cargar-datos

    set pob-prev 0


    set pob-prev count patches with [pcolor != white and pxcor = 0]


    if (activar-penalizacion = True) [penalizar-similares]
    ejecutar-regla-EIA-AB

    ;show pher-prom

    ifelse (ticks > 100) [
      set pher-prom fput pheromones pher-prom
      set pher-prom but-last pher-prom

    ]
    [
      set pher-prom fput pheromones pher-prom
    ]



    ;show pher-prom
    ;show "--"


    if(ticks > 100) [strong-agent]
    if(ticks > 25) [calc-mean-energy]
    if (ticks > 300) [ inoculate-agent ]
    if (ticks > 500) [check-attack]
    delta-poblacion
    registro-funcionamiento
    calculo-umbral
    graficar
    verificar-ataque
    tabla-contingencia


    ;;HISTOGRAMA;;;;;;;;;;;;;;;;;;;;;;;;
    calcular-diversidad               ;;
    preparar-histograma               ;;
    graficar-diversidad               ;;
    do-plots4                         ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    tick
  ]
  [
    stop
  ]
end

to strong-agent
  ;show "--"
 ; show mean pher-prom
  if (mean pher-prom > 0 and mean pher-prom < 0.15)[
    ;show "Getting strong agent"
    set flag_capsule 1
     ask max-n-of 50 patches with [pxcor = 0] [energia][
      set cryogenic-capsule lput genes cryogenic-capsule
    ]
  ]


end

to calc-mean-energy

  if (mean energy_prom > energy_umbral)[set energy_umbral (mean energy_prom)]

end

to check-attack

  ifelse ((energy_umbral - mean (energy_prom)) > energy_epsylon ) [
    set timer_attack (timer_attack + 1)
    if(timer_attack > timer_umbral) [
      ;show "UNDER ATTACK!"
      set under_attack 1
      set timer_normality 0
    ]
  ]
  [
    set timer_normality (timer_normality + 1)
    if(timer_normality > 25)[
      ;show "NORMALITY"
      set under_attack 0
      set timer_attack 0
    ]
  ]


end



to inoculate-agent

  if(pob-delta < -5 and flag_capsule = 1)[
    show "Inoculate"
    foreach cryogenic-capsule[ i ->
      ask one-of ( patches with [pxcor = 0] ) [
        set genes i
      ]
    ]
  ]



end

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Rutinas de Visualización de comportamiento
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



  to do-plots4
  set-current-plot "Variabilidad Genética" ;; which plot we want to use next
  set-current-plot-pen "diversidad" ;; which pen we want to use next
  plot  indicador-variabilidad ;; what will be plotted by the current pen

end






  to graficar

  set celdas-verdes count patches  with [pxcor = 0 and pcolor = green]
  set celdas-amarillas count patches  with [pxcor = 0 and pcolor = yellow]
  set celdas-rojas count patches  with [pxcor = 0 and pcolor = red]


  set energia-acumulado 0
  ask patches with [pxcor = 0] [set energia-acumulado (energia-acumulado + energia) ]
  set energia-promedio-algas (energia-acumulado / PAR-cantidad-celdas)

  set-current-plot "Celdas-Tiempo"
  ;;set-current-plot-pen "celdas-verdes"
  ;;plot celdas-verdes
  ;;set-current-plot-pen "celdas-amarillas"
  ;;plot celdas-amarillas
  ;;set-current-plot-pen "celdas-rojas"
  ;;plot celdas-rojas
  set celdas-totales (celdas-verdes + celdas-amarillas + celdas-rojas )
  set-current-plot-pen "celdas-totales"
  plot celdas-totales
  set-current-plot-pen  "energia-promedio"
  plot energia-promedio-algas
  ;;set-current-plot-pen  "banda-superior-celdas"
  set-current-plot-pen  "banda-inferior-celdas"
  plot banda-inferior-poblacion
  ;;set-current-plot-pen  "banda-superior-energia"
  set-current-plot-pen  "banda-inferior-energia"
  plot banda-inferior-energia
  set-current-plot "Alertas-Tiempo"
  set-current-plot-pen "umbral"
  plot umbral-alerta
  set-current-plot-pen "alertas"
  plot salida-banda



  end

  to calcular-diversidad

  set variable-histograma []

  ask patches with [pxcor = 0 and energia > 0]  [
              let i 0
              let caracteristica  0
              let representante 1
              let caracteristica-primo 0


              while[i < 5] ;;no puede ser más grande seguramente por problemas de precisión matemática de herramienta.
                   [
                   ;;deteminar el gen de la posici�n i
                   set caracteristica item i genes
                   set caracteristica-primo item caracteristica numeros-primos

                   ;;calcular representante
                   set representante (representante * (caracteristica-primo ^ i ) )
                   set i (i + 1)

                  ]
                  set representante round (abs (representante mod 2389 ))
                  ;; set representante (1 )
                  set variable-histograma lput representante variable-histograma

            ]

end




  to graficar-diversidad

      histogram variable-histograma

;;------------------------CALCULO DE INDICE DE VARIABILIDAD--------------------------

let temporal-histograma variable-histograma
let pre-descarte 0
let post-descarte 0
set numero-genotipos 0
set poblaciones []
let suma-inferior 0


while[ (not empty? temporal-histograma )]
     [
     ;; conteo pre-eliminaci�n de la lista
     set pre-descarte length temporal-histograma
     let a (item 0  temporal-histograma) set temporal-histograma remove a temporal-histograma set numero-genotipos (numero-genotipos + 1) set flag-diversidad true
     ;; conteo post-eliminacion de la lista
     set post-descarte length temporal-histograma


     ;;agregar a lista
     set poblaciones lput ( pre-descarte - post-descarte ) poblaciones
     ]

    ;;realizar calculo de indice de diversidad y entregar

    ;;sumando acumulados
    while [ (not empty? poblaciones) ]

    [
    set suma-inferior (suma-inferior + first poblaciones)
    set poblaciones remove-item 0 poblaciones
    ]
    if (suma-inferior > 0 ) [set indicador-variabilidad  numero-genotipos / suma-inferior]
end

to preparar-histograma
  set-current-plot "Variabilidad-Genetica"
  set-plot-x-range 0 2389
  set-plot-y-range 0 (100)
  set-histogram-num-bars longitud-automata


end

to registro-funcionamiento

  let i max-pxcor
  repeat (max-pxcor) [
    ask patches with [pxcor = i] [
      set pcolor [pcolor] of patch-at -1 0
    ]
    set i (i - 1)
  ]



end



to ejecutar-regla-EIA-AB
  ; TO-DO
  ; si tiene demasiada ventaja, mata al vecino y pone su clon

  ask patches with [pxcor = 0] [
    let genP genes
    let aux 0
    let aux2 0
    let rand 0
    let rand2 0
    set energia (energia - 1) ; consumo energia basal
    if (energia < 0) [ set energia 0 ]

    if (energia > energia-maxima) [

      let cont_rep 1

      if ([energia] of patch-at 0 1 = 0)[
        ask patch-at 0 1 [
          set energia (energia-maxima * 0.5)
          set rand random 53
          set rand2 random 53
          set aux item rand genP
          set aux2 item rand2 genP
          set genP replace-item rand genP aux2
          set genP replace-item rand2 genp aux
          set genes genP
          set cont_rep cont_rep + 1
        ]
      ]
      if ([energia] of patch-at 0 -1 = 0)[
        ask patch-at 0 -1 [
          set energia (energia-maxima * 0.5)
          set rand random 53
          set rand2 random 53
          set aux item rand genP
          set aux2 item rand2 genP
          set genP replace-item rand genP aux2
          set genP replace-item rand2 genp aux
          set genes genP
          set cont_rep cont_rep + 1
        ]
      ]
      ;if ([energia] of patch-at 0 1 > 0)[ask patch-at 0 1 [set energia energia-maxima]]
      ;if ([energia] of patch-at 0 1 > 0)[ask patch-at 0 1 [set energia (energia + ([energia] of myself) * 0.5)]]

      set energia energia-maxima / 3
    ]



    if (energia > 100) [ set pcolor green ]
    if ((energia <= 100) and ( energia > 50)) [ set pcolor yellow ]
    if ((energia <= 50) and ( energia > 0)) [ set pcolor red ]
    if (energia <= 0) [
      set pcolor white
      set pheromones pheromones + 1
    ]
    if (pheromones > 0)[
      set pheromones  pheromones * 6 / 7
      if (pheromones < 0) [
        set pheromones 0
      ]
    ]

  ]
end

to penalizar-similares

  let AtGen n-values 53 [0]
  let LowGen 0
  let Comparation 0
  let Matches 0
  ask patches with [ pxcor = 0][
    set AtGen genes
    set LowGen [genes] of patch-at 0 -1
    set Comparation (map > AtGen LowGen)
    set Matches occurrences false Comparation

    ;show Matches

    if (Matches > 50)[
      ;show "energia anterior:"
      ;show energia
      set energia (energia * 0.5)
      ;show "energia actual:"
      ;show energia
    ]
  ]



  ; TO-DO
end

to-report occurrences [x the-list]
  report reduce
    [ [occurrence-count next-item] -> ifelse-value (next-item = x) [occurrence-count + 1] [occurrence-count] ] (fput 0 the-list)
end


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Rutina de carga de datos desde archivo
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to cargar-datos


  if(Tipo_ataque = "variedLength")[
    ;show "Distintos Ataques"

  ]

  if(Tipo_ataque = "superSet")[
    ;show "Distintos Ataques"
    superSet
  ]

  if(Tipo_ataque = "Distintos Ataques")[
    ;show "Distintos Ataques"
    cargar-datos2
  ]


  if(Tipo_ataque = "Ataques Esporadicos")[
    ;show "Ataques Esporadicos"
    cargar-datos3
  ]




end


to variedLength





end


to superSet

  if ( (flag-par = 1) and (not (file-at-end?))  ) [ ;; hay lineas en el archivo que leer


    set lineas-leidas lineas-leidas + 1

    ;; se trabaja como string se deben individualizar los datos y usar reglas para generar breeds sedimientos (en las reglas se puede manejar estados multiples y otras discretizaciones
    ;; luego de generar los sedimentos "nuevos" se procede con las herramientas de an�lisis
    ;; eventualmente cambiar preNN.pl para facilitar la lectura de datos



    set caracteristicas-archivo  (sentence    list file-read file-read  file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read
      file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read
      file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read)

    ifelse ((item 54 caracteristicas-archivo) > 0) [set ataque 400 ] [set ataque 0]

    ;; Esta versión carga las características que son binarias al flujo
    ;;"puertoOrigen_0,puertoDestino_1,protocolo_2,TTL_3,TOS_4,IPLen_5,DgmLen_6,RB_7,MF_8,DF_9,opcionesIP_10,F1_11,F2_12,U_13,A_14,P_15,R_16,S_17,F_18,Win_19,
    ;;TcpLen_20,opcionesTCP_21,UDPLen_22,Type_23,Code_24,telnet_25,ssh_26,ftp_27,netbios_28,rlogin_29,rpc_30,nfs_31,lockd_32,netbiosWinNT_33,Xwin_34,dns_35,
    ;;ldap_36,smtp_37,pop_38,imap_39,http_40,ssl_41,px_42,serv_43,time_44,tftp_45,finger_46,nntp_47,ntp_48,lpd_49,syslog_50,snmp_51,bgp_52,socks_53\n"

    ;; configuraci�n binaria
    set indice-caracteristicas-utilizadas [7 8 9 11 12 13 14 15 16 17 18 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53]


    ;; determinando el número de características presentes (binario)

    let contador-caracteristicas-activas 0
    set particulas[]
    foreach indice-caracteristicas-utilizadas [ ?1 ->

      if ((item ?1 caracteristicas-archivo) > 0)
      [
        set contador-caracteristicas-activas (contador-caracteristicas-activas + 1)
        set particulas lput ?1 particulas



      ]

    ]

    ;;Alimentar celdas (agentes)

    ask patches with [ pxcor = 0 and pcolor != white][
      foreach particulas [ i ->
        set energia (energia + 80 - ((position i genes) ^ 2) * 0.092 )
        ; 0.087
        if ([energia] of patch-at 0 1 > 0)[ask patch-at 0 1 [set energia (energia + (([energia] of myself) * 0.1))]]
        if ([energia] of patch-at 0 -1 > 0)[ask patch-at 0 -1 [set energia (energia + (([energia] of myself) * 0.1))]]
      ]

    ]


      ifelse (ticks > 25) [
        set energy_prom fput energia-promedio-algas energy_prom
        set energy_prom but-last energy_prom

      ]
      [
        set energy_prom fput energia-promedio-algas energy_prom
      ]


    ;; el orden de los datos importa mucho

    ;;Aplicación del metabolismo

  ] ;;cierre de if de verificación de archivo

  if ((flag-par = 1) and index_attack > 19)
  [
      set fin-simulacion true
      stop
  ]

  if((flag-par = 1) and ( ticks > (  item index_attack List_ss + 2000)  ))[
    show "Back to normal!"
    set index_attack index_attack + 1
    set fuente-datos-bin 0
  ]
  if((flag-par = 1) and ( ticks =  (  item index_attack List_ss )))[
  show "Under Attack!"
  set fuente-datos-bin 1
  ]

  if((flag-par = 1) and ticks = 3905)[
    show "End of Training"
  ]

    ;;Primera tanda de ataques de Denegacion de Servicio
    if ( flag-par = 0 )
    [
      set nombre-archivo archivoAtaque

      file-open nombre-archivo
      set fuente-datos nombre-archivo
      set fuente-datos-bin 0
      set archivos-procesados (archivos-procesados + 1)
      set lineas-leidas 0
      set flag-par 1

    ]



end


to cargar-datos2  ;;version con carga de datos desde archivo
  ifelse ( (lineas-leidas < 950) and (not (file-at-end?))  ) [ ;; hay lineas en el archivo que leer


    set lineas-leidas lineas-leidas + 1

    ;; se trabaja como string se deben individualizar los datos y usar reglas para generar breeds sedimientos (en las reglas se puede manejar estados multiples y otras discretizaciones
    ;; luego de generar los sedimentos "nuevos" se procede con las herramientas de an�lisis
    ;; eventualmente cambiar preNN.pl para facilitar la lectura de datos



    set caracteristicas-archivo  (sentence    list file-read file-read  file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read
      file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read
      file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read)

    ifelse ((item 54 caracteristicas-archivo) > 0) [set ataque 400 ] [set ataque 0]

    ;; Esta versión carga las características que son binarias al flujo
    ;;"puertoOrigen_0,puertoDestino_1,protocolo_2,TTL_3,TOS_4,IPLen_5,DgmLen_6,RB_7,MF_8,DF_9,opcionesIP_10,F1_11,F2_12,U_13,A_14,P_15,R_16,S_17,F_18,Win_19,
    ;;TcpLen_20,opcionesTCP_21,UDPLen_22,Type_23,Code_24,telnet_25,ssh_26,ftp_27,netbios_28,rlogin_29,rpc_30,nfs_31,lockd_32,netbiosWinNT_33,Xwin_34,dns_35,
    ;;ldap_36,smtp_37,pop_38,imap_39,http_40,ssl_41,px_42,serv_43,time_44,tftp_45,finger_46,nntp_47,ntp_48,lpd_49,syslog_50,snmp_51,bgp_52,socks_53\n"

    ;; configuraci�n binaria
    set indice-caracteristicas-utilizadas [7 8 9 11 12 13 14 15 16 17 18 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53]


    ;; determinando el número de características presentes (binario)

    let contador-caracteristicas-activas 0
    set particulas[]
    foreach indice-caracteristicas-utilizadas [ ?1 ->

      if ((item ?1 caracteristicas-archivo) > 0)
      [
        set contador-caracteristicas-activas (contador-caracteristicas-activas + 1)
        set particulas lput ?1 particulas



      ]

    ]

    ;;Alimentar celdas (agentes)

    ask patches with [ pxcor = 0 and pcolor != white][
      foreach particulas [ i ->
        set energia (energia + 83 - position i genes * 2.9)
      ]


    ]
    ;; el orden de los datos importa mucho

    ;;Aplicación del metabolismo

  ] ;;cierre de if de verificación de archivo

  [
    if (flag-par = 1)
    [
      set rounds (rounds + 1)
      ifelse rounds >= 2 [
        set fin-simulacion true
        stop
      ]
      [
        file-open "datosNormalesEIA3000.dat"
        set fuente-datos "datosNormalesEIA3000.dat"
        set nombre-archivo "datosNormalesEIA3000.dat"
        set lineas-leidas 0
        ifelse(Tipo_ataque_value = 0)
      [
        set Tipo_ataque_value 1
      ]
      [
        set Tipo_ataque_value 0
      ]
        show "Back to normal!"
      ]
    ]

    ;;Primera tanda de ataques de Denegacion de Servicio
    if ( flag-par = 0 )
    [
      ifelse(Tipo_ataque_value = 0)
      [
        set nombre-archivo archivoAtaque
      ]
      [
        set nombre-archivo archivoAtaque2
      ]
      show nombre-archivo

      file-open nombre-archivo
      set fuente-datos nombre-archivo
      set fuente-datos-bin 1
      set archivos-procesados (archivos-procesados + 1)
      set lineas-leidas 0

    ]

    ifelse (flag-par = 0 )
    [set flag-par 1]
    [set flag-par 0]
  ]
end

to cargar-datos3  ;;version con carga de datos desde archivo -- Ataques Esporadicos
  ifelse ( (lineas-leidas < 300) and (not (file-at-end?))  ) [ ;; hay lineas en el archivo que leer


    set lineas-leidas lineas-leidas + 1

    ;; se trabaja como string se deben individualizar los datos y usar reglas para generar breeds sedimientos (en las reglas se puede manejar estados multiples y otras discretizaciones
    ;; luego de generar los sedimentos "nuevos" se procede con las herramientas de an�lisis
    ;; eventualmente cambiar preNN.pl para facilitar la lectura de datos



    set caracteristicas-archivo  (sentence    list file-read file-read  file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read
      file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read
      file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read file-read)

    ifelse ((item 54 caracteristicas-archivo) > 0) [set ataque 400 ] [set ataque 0]

    ;; Esta versión carga las características que son binarias al flujo
    ;;"puertoOrigen_0,puertoDestino_1,protocolo_2,TTL_3,TOS_4,IPLen_5,DgmLen_6,RB_7,MF_8,DF_9,opcionesIP_10,F1_11,F2_12,U_13,A_14,P_15,R_16,S_17,F_18,Win_19,
    ;;TcpLen_20,opcionesTCP_21,UDPLen_22,Type_23,Code_24,telnet_25,ssh_26,ftp_27,netbios_28,rlogin_29,rpc_30,nfs_31,lockd_32,netbiosWinNT_33,Xwin_34,dns_35,
    ;;ldap_36,smtp_37,pop_38,imap_39,http_40,ssl_41,px_42,serv_43,time_44,tftp_45,finger_46,nntp_47,ntp_48,lpd_49,syslog_50,snmp_51,bgp_52,socks_53\n"

    ;; configuraci�n binaria
    set indice-caracteristicas-utilizadas [7 8 9 11 12 13 14 15 16 17 18 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53]


    ;; determinando el número de características presentes (binario)

    let contador-caracteristicas-activas 0
    set particulas[]
    foreach indice-caracteristicas-utilizadas [ ?1 ->

      if ((item ?1 caracteristicas-archivo) > 0)
      [
        set contador-caracteristicas-activas (contador-caracteristicas-activas + 1)
        set particulas lput ?1 particulas



      ]

    ]

    ;;Alimentar celdas (agentes)

    ask patches with [ pxcor = 0 and pcolor != white][
      foreach particulas [ i ->
        set energia (energia + 83 - position i genes * 2.9)
      ]


    ]
    ;; el orden de los datos importa mucho

    ;;Aplicación del metabolismo

  ] ;;cierre de if de verificación de archivo

  [
    if (flag-par = 1)
    [
      set rounds (rounds + 1)
      ifelse rounds >= 6 [
        set fin-simulacion true
        stop
      ]
      [
        file-open "datosNormalesEIA3000.dat"
        set fuente-datos "datosNormalesEIA3000.dat"
        set nombre-archivo "datosNormalesEIA3000.dat"
        set lineas-leidas 0
        show "Back to normal!"
      ]
    ]

    ;;Primera tanda de ataques de Denegacion de Servicio
    if ( flag-par = 0 )
    [
      set nombre-archivo archivoAtaque

      file-open nombre-archivo
      set fuente-datos nombre-archivo
      set fuente-datos-bin 1
      set archivos-procesados (archivos-procesados + 1)
      set lineas-leidas 0

    ]

    ifelse (flag-par = 0 )
    [set flag-par 1]
    [set flag-par 0]
  ]
end





  to calculo-umbral
    ;;ingresando nueva energía en la lista

    ;;-------------------------------------Implementar Energía--------------------------------------
    set variable-ventana-deslizante-energia lput energia-promedio-algas variable-ventana-deslizante-energia
    set variable-ventana-deslizante-energia but-first variable-ventana-deslizante-energia



    set variable-ventana-deslizante-poblacion lput celdas-totales variable-ventana-deslizante-poblacion
    set variable-ventana-deslizante-poblacion but-first variable-ventana-deslizante-poblacion
    ;;-------------------------------------Cálculo de media y desviación estándar-------------------
    ;;ENERGIA
    set media-muestral-deslizante-energia mean variable-ventana-deslizante-energia
    set desviacion-muestral-deslizante-energia standard-deviation  variable-ventana-deslizante-energia
    set banda-superior-energia (media-muestral-deslizante-energia + PAR-banda-superior-energia )
    set banda-inferior-energia (media-muestral-deslizante-energia - PAR-banda-inferior-energia )

    ;;POBLACION
    set media-muestral-deslizante-poblacion mean variable-ventana-deslizante-poblacion
    set desviacion-muestral-deslizante-poblacion standard-deviation  variable-ventana-deslizante-poblacion
    set banda-superior-poblacion (media-muestral-deslizante-poblacion + PAR-banda-superior-poblacion )
    set banda-inferior-poblacion (media-muestral-deslizante-poblacion - PAR-banda-inferior-poblacion )

    ;;consideración de límites

    if (banda-inferior-poblacion < limite-poblacion )[ set banda-inferior-poblacion limite-poblacion]
    if (banda-inferior-energia < limite-energia )[ set banda-inferior-energia limite-energia]



  end

to verificar-ataque

  ifelse((energia-promedio-algas < banda-inferior-energia) or (celdas-totales <

    banda-inferior-poblacion))
  [set salida-banda salida-banda + factor-aumento-alerta]
  [if (salida-banda > 0 ) [set salida-banda salida-banda - factor-disminucion-alerta]]


  ifelse (salida-banda > umbral-alerta )
  [set ataque-detectado true
    if ticks > 950 and ticks < 1900 and activeDelay1 = 0[
      set activeDelay1 (ticks - 950)
    ]
    if ticks > 2850 and activeDelay2 = 0[
      set activeDelay2 (ticks - 2850)
    ]
  ]
  [set ataque-detectado false]


end



to tabla-contingencia
  if ( ticks > ventana-entrenamiento)
  [
    ;;if (ataque-detectado) [set contador-alertas (contador-alertas+1)]
    set entrenamiento "clasificando"
    ;;Momento para intervenir si se considera congelar la evolución

    ;;Muestreo Unitario, válido para calibración de parámetros basado en paquete a paquete
    if (( fuente-datos-bin = 1 ) and (ataque-detectado)) [set VP (VP + 1)]
    if (( fuente-datos-bin = 1 ) and (not ataque-detectado)) [set FN (FN + 1)]
    if (( fuente-datos-bin = 0 ) and (not ataque-detectado)) [set VN (VN + 1)]
    if (( fuente-datos-bin = 0 ) and (ataque-detectado)) [set FP (FP + 1)]


    ;; Calculando Sensibilidad y Especificidad

    if ( (VP + FN) > 0 )
    [set sensibilidad  (VP / (VP + FN) )]
    if ( (FP + VN) > 0 )
    [set especificidad (VN / (FP + VN) )]




  ]




end

to establecer-limites-bandas
  let media-limite-energia mean variable-ventana-deslizante-energia
  let desviacion-limite-energia standard-deviation  variable-ventana-deslizante-energia
  let media-limite-poblacion mean variable-ventana-deslizante-poblacion
  let desviacion-limite-poblacion standard-deviation  variable-ventana-deslizante-poblacion

  set limite-energia  (media-limite-energia -  desviacion-limite-energia)
  set limite-poblacion (media-limite-poblacion -  desviacion-limite-poblacion)


end
@#$#@#$#@
GRAPHICS-WINDOW
711
442
1321
653
-1
-1
2.0
1
10
1
1
1
0
0
1
1
0
300
0
100
1
1
1
ticks
30.0

TEXTBOX
725
422
875
440
AC
11
0.0
1

BUTTON
577
26
659
59
Inicializar
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
577
66
655
99
Ejecutar
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
657
181
1320
435
Celdas-Tiempo
Tiempo
#
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"celdas-totales" 1.0 0 -14439633 true "" ""
"energia-promedio" 1.0 0 -3844592 true "" ""
"band-superior-celdas" 1.0 0 -11221820 true "" ""
"banda-inferior-celdas" 1.0 0 -13791810 true "" ""
"banda-superior-energia" 1.0 0 -1264960 true "" ""
"banda-inferior-energia" 1.0 0 -3508570 true "" ""
"umbral energia" 1.0 0 -1184463 true "" "plot energy_umbral"
"prom-energia-actual" 1.0 0 -16777216 true "" "plot mean energy_prom"

PLOT
3
12
232
253
Variabilidad-Genetica
Especímenes
Cantidad
0.0
10.0
-1.0
200.0
true
false
"" ""
PENS
"default" 1.0 1 -5825686 true "" ""

PLOT
3
336
494
486
Variabilidad Genética
Tiempo
Indicador de Variabilidad
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"diversidad" 1.0 0 -817084 true "" ""

SWITCH
53
580
218
613
activar-penalizacion
activar-penalizacion
1
1
-1000

MONITOR
524
234
648
279
Muestras Revisadas
muestras-revisadas
17
1
11

MONITOR
502
339
647
384
Situación
fuente-datos
17
1
11

MONITOR
524
284
624
329
Entrenamiento?
Entrenamiento
17
1
11

TEXTBOX
497
106
647
124
Tabla de Contingencia
11
0.0
1

MONITOR
528
136
585
181
NIL
VP
17
1
11

MONITOR
588
136
645
181
NIL
FP
17
1
11

MONITOR
528
183
585
228
NIL
FN
17
1
11

MONITOR
589
183
646
228
NIL
VN
17
1
11

PLOT
711
682
1322
832
Alertas-Tiempo
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"alertas" 1.0 0 -10022847 true "" ""
"umbral" 1.0 0 -16448764 true "" ""

MONITOR
498
397
573
442
NIL
sensibilidad
17
1
11

MONITOR
500
455
583
500
NIL
especificidad
17
1
11

MONITOR
592
474
704
519
NIL
fuente-datos-bin
17
1
11

MONITOR
590
524
677
569
NIL
limite-energia
17
1
11

MONITOR
488
525
585
570
NIL
limite-poblacion
17
1
11

BUTTON
508
66
572
99
step
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
241
512
484
572
archivoAtaque
superSet01
1
0
String

PLOT
788
10
1268
179
Pheromones
NIL
NIL
0.0
0.0
0.0
0.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot pheromones"

INPUTBOX
241
585
484
645
ArchivoAtaque2
Tor-00
1
0
String

CHOOSER
60
517
233
562
Tipo_ataque
Tipo_ataque
"variedLength" "superSet" "Ataques Esporadicos" "Distintos Ataques"
1

PLOT
274
72
474
222
Variedad energia
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count turtles"

INPUTBOX
655
113
737
173
energy_epsylon
20.0
1
0
Number

INPUTBOX
660
46
768
106
timer_umbral
100.0
1
0
Number

MONITOR
461
235
518
280
attack
under_attack
17
1
11

@#$#@#$#@
#EL MODELO
Esta versión del modelo, tal como describe el proyecto, pretende ser una simplificación del modelo EIA original apuntando a lograr un algoritmo más limpio y eficiente que funcione como clasificador. (leer el proyecto)

##CARACTERÍSTICAS PRINCIPALES DEL MODELO

Este modelo migra desde los ABM (Agents Based Models) al terreno de los CA (Cellular Automata), pues en este caso el modelo no tiene un conjunto de agentes que se reproducen en una topología bidimensional, sino que  un conjunto de celdas en un anillo (autómata celular lineal), por lo cual la cantidad de agentes es fijo. Desde este punto de vista, no se crean nuevos agentes, sino que las celdas pueden cambiar el estado de sus vecinas  cuando su energía supera al de estas por una tasa dada, manteniendo para esto las reglas de mutación del modelo original.

Todas las celdas reciben todas las partículas, y no aleatoreamente con en la versión anterior.

Se sigue midiendo la energía promedio y la cantidad de celdas vivas ( celdas del autómata con energía).
Se agregó una heurísticas de penalización de vecinos similares para mejorar la variabilidad genética del modelo, con lo cual se busca que el modelo se ajuste a todo el escenario de "normalidad" y no sólo a condiciones triviales.

###Otros Cambios
En la determinación de alertas mediante umbral ya no se realiza una acumulación lineal de instancias de salida de bandas hasta traspasar el umbral máximo para generar una alerta y luego volverse cero nuevamente con el retorno a la banda del estado de energía o tamaño de población, sino que existe parámetros de  "incremento" y "decremento" para la salida y el retorno de bandas.

Los parámetros han sido calibrados usando Metaaprendizaje con el uso de la herramienta Behavior Search (http://behaviorsearch.org/) mediante un algoritmo genético.


###Sobre la Interfaz

La interfaz principal del modelo presenta un rectángulo (AC) donde en la posición más extrema de la izquierda se representa el estado actual del AC de celdas (vertical, cada fila es una celda) y hacia la derecha la historia del mismo en el tiempo.

El modelo provee un histograma de "biotipos" de celdas y un gráfico que contabiliza la variabilidad genética total, en el cual se puede apreciar como decae en el tiempo.


#EL PROYECTO

Formalización y ajuste de parámetros del Modelo de Especies Indicadoras Artificiales (EIA), como Clasificador aplicado a la Detección de Intrusos en Redes.

## RESUMEN

El objetivo del proyecto consiste en mejorar el modelo de Especies Indicadoras Artificiales (EIA), logrando con esto resultados competitivos con  respecto a algoritmos clasificadores en el estado del arte de la detección de intrusos en redes de computadoras.
El modelo EIA [1,2], corresponde al desarrollo de un  Sistema de Detección de Anomalías en flujo (SDAF), basado en  conceptos a Vida Artificial [3] y Sistemas Inmunes Artificiales [4], que ha demostrado  tener buenos resultados comparados con otras técnicas de machine learning en términos de reducción de falsos positivos y de exactitud, para escenarios de prueba basados en Benchmark [5]  estándar del área.
En términos simples el modelo EIA corresponde al desarrollo de una algoritmo clasificador que aprende de forma autónoma a diferenciar un flujo de eventos normales de los que son anormales.

Este trabajo se hace necesario debido a que el actual modelo EIA ha sido desarrollado en base a herramientas de prototipado de alto nivel, no teniendo la eficiencia necesaria para ser puesto a nivel de producción. Por lo cual,  este proyecto busca desarrollar un algoritmo computacionalmente eficiente en base al modelo EIA.

El proyecto de mejoramiento del modelo EIA contempla la simplificación de este tomando como base una migración del dominio de Modelos Basados en Agentes (ABM)  al de Autómatas Celulares (CA). Por otro lado, una vez migrado el modelo, se considera la aplicación de técnicas de Meta-Aprendizaje para la calibración de parámetros del modelo.

Las Técnica de meta-aprendizaje están referidas a sistemas híbridos donde, más de una metaheurísticas es probada usualmente en un diseño de capas. En este caso  un algoritmo genético tutor será usado para buscar la mejor configuración de parámetros para el modelo.

Esta investigación será realizada en colaboración con investigadores de la Universidad Politécnica de Madrid.

Como resultado de esta investigación se espera ponencias para dos conferencias internacionales y un artículo para revista indexada (ISI). 

Referencias:
[1] Pinacho Pedro, Especies Indicadoras Artificiales (EIA) Una técnicas de detección de anomalías en flujo basada en vida artificial, Tesis de Magíster en Ingeniería Informática, Departamento de Ingeniería Informática, Universidad de Santiago de Chile, 2011.
[2] Pinacho Pedro, Pau Iván, Chacón Max y Sánchez Sergio, An ecological approach to anomaly detection: the EIA model, ICARIS’12, Taormina, Italy, August 2012.
[3] Langton, C.G., Artificial life : an overview. Complex adaptive systems. 1995, Cambridge, Mass.: MIT Press. xi, 340 p., [6] p. of plates.
[4] Glickman, M., J. Balthrop, and S. Forrest, A machine learning evaluation of an artificial immune system. Evolutionary Computation, 2005. 13(2): p. 179-212.
[5] MIT. DARPA Intrusion Detection Evaluation.  1998; Datos para Benchmarking de sistemas IDS. 


##OBJETIVOS
###Objetivo General
Perfeccionar el modelo Especies Indicadoras Artificiales (EIA), logrando  un algoritmo  computacionalmente eficiente que demuestre ser  superior a otros clasificadores en el estado del arte en la detección de ataques innovadores.

###Objetivos Específicos
1-Estudiar la aplicación de nuevas heurísticas biológicas al modelo actual, y verificar su impacto.
2-Simplificar y formalizar  el Modelo EIA usando el paradigma de Autómatas Celulares.
3-Realizar un ajuste de parámetros usando técnicas de meta-aprendizaje por medio de un algoritmo genético.
4-Diseñar una plan de pruebas para el nuevo clasificador y otros algoritmos en el estado del arte, para comparación.
5-Desarrollar pruebas comparativas y capturar resultados
6-Verificar potencial reactivo del algoritmo (*)
7-Análisis de resultados y escribir reportes.

(*) Objetivo secundario, pensado en lograr que el modelo no sólo opere como clasificador sino que logre intervenir su entorno para mantener la estabilidad, esto es actué como  un Sistema de Prevención de Intrusiones (IPS). O IDS Reactivo.


##HIPÓTESIS

Hipótesis: El modelo EIA puede mejorar su bondad como sistema clasificador y eficiencia computacional a través de un proceso de simplificación y formalización matemática, logrando superar a clasificadores en el estado del arte en el problema de detección de intrusos en redes.

##METODOLOGÍA

#Metodología

La metodología de desarrollo, contempla el desarrollo de prototipos desechables y prototipos incrementales para apoyar el método científico [16].

###Datos de Prueba
Las pruebas del EIA son realizadas a través del Set de Datos DARPA’ 98, descartándose el uso de los Datos KDD’99, por ser considerados perjudiciales. KDD’99 fue creado a partir de una porción de datos del set DARPA’98 ampliamente discretizada, perdiendo con esto validez y utilidad. Siendo no recomendado su uso ni aceptación para la validación de técnicas de detección de intrusos [11, 61].
DARPA’98
El grupo de Tecnología de Sistemas de Información (IST) del Lincoln Laboratory del MIT durante los años 1998 y 1999, bajo el amparo de la agencia de proyectos avanzados de investigación de defensa (DARPA ITO) y el patrocinio del laboratorio de investigación de la Fuerza Aérea de EEUU, recolectó y distribuyó el primer conjunto estándar de datos para la evaluación de IDS en redes computacionales.
Este conjunto de datos corresponde a una simulación de una red de área local (LAN, Local Area Network) de una base aérea estadounidense, en la cual ocurrían ataques a diferentes estaciones de trabajo. El resultado, fueron semanas de datos en bruto de las conexiones de red y el etiquetado de las sesiones con los diferentes tipos de ataques.
El conjunto de datos del año 1998, se compone de  nueve semanas de datos, a su vez cada semana por cinco días (Lunes a Viernes). Las nueve semanas se dividen en siete semanas de datos para entrenamiento y dos semanas de datos para prueba.

Este set de datos provee decenas de miles de paquetes de red para análisis, los cuales serán utilizados para las pruebas consideradas.

Las comparaciones de los algoritmos serán realizadas en los términos habituales para este tipo de algoritmos: tablas de contingencia (falsos positivos (FP) , falsos negativos (FN), verdaderos positivos (VP), verdaderos negativos (VN)), además de sensibilidad, especificidad y exactitud.

###Validación de significancia estadística

La validación de la significancia se realizará a través del Test de los signos de Wilcoxon.
El test de los signos de Wilcoxon es una prueba no paramétrica para comparar la mediana de dos muestras relacionadas y determinar si existen diferencias entre ellas. Se utiliza como alternativa a la prueba  t de Student cuando no se puede suponer la normalidad de dichas muestras. Debe su nombre a Frank Wilcoxon, que la publicó en 1945.
Se utiliza cuando la variable subyacente es continua pero no presupone ningún tipo de distribución particular.
Con el presente test se comparan los clasificadores en término de la exactitud y  tasa de falsos positivos. Si el valor de p es mayor a 0,05 quiere decir que los clasificadores comparados no tienen diferencias significativas. Al contrario, si p es menor a 0,05 quiere decir que tienen diferencias significativas.






## FUNDAMENTACIÓN DEL PROYECTO - MARCO TEÓRICO

###El problema

El problema abordado en este proyecto está enmarcado dentro del área de los  Sistema de detección de intrusos (IDS).  Estos sistemas son herramientas tecnológicas que  aportan una capacidad de detección al conjunto de defensas informáticas de un sistema, alertando acerca de la presencia de actividad sospechosa que ocurra antes, durante  o posterior a un ataque [7].
Existen diversas formas de clasificar  estos sistemas, dependiendo de sus características, determinando en base a esto el tipo particular de problema a abordar [3], entre estas tenemos:

•Técnica de  Detección: Un sistema IDS puede usar técnicas de detección basada en firmas o patrones (conocimiento experto), o bien utilizar alguna técnica para detección de anomalías. Este trabajo se centra en esta segunda forma.
•Fuente de Datos: Para la realización de la clasificación el sistema IDS debe analizar datos en flujo, llamados también muchas veces audit trails, los cuales pueden ser entendidos como una secuencia de registros de eventos ordenados temporalmente que manifiestan la actividad de usuarios o recursos del sistema. Tradicionalmente [8] los sistemas pueden tomar sus datos de registros de eventos sintetizados por los sistemas computacionales o Hosts, por lo cual reciben el nombre de HIDS o de la red de forma directa a través de técnicas de sniffing, entendiendo esto como una forma de escuchar la red de forma promiscua (todo el tráfico sin importar el destinatario), estos sistemas se suelen llamar NIDS (Network Intrusion Detection Systems), estando este trabajo enmarcado en este tipo de propuestas.
•Reactividad: Los sistemas de detección de intrusos pueden actuar de formar reactiva, no sólo informando la presencia de una intrusión detectada sino que además tomando medidas para tratar de neutralizarla. Este es el caso de los actuales sistemas IPS (Intrusion Prevention System), los cuales de igual forma que sistema Firewall pueden tomar medidas como  por ejemplo el cierre de sesiones hostiles. Debido a que este trabajo se orienta a la obtención de una nueva técnica de clasificación de datos en flujo y en particular al análisis de tráfico de red para la detección o clasificación de amenazas, la reactividad está fuera del área de estudio, acotando el proyecto a la visión tradicional de IDS, como ente de detección y reporte de ataques; no obstante esto, se considera en este proyecto el estudio del protencial reactivo del modelo EIA.

##Problemas de las soluciones actuales

Tanto las propuestas de detección de intrusos basadas en conocimiento experto como las basadas en técnicas de aprendizaje automático tradicional tales como Redes Neuronales [9], Support Vector Machine [5] entre otras poseen una desventaja operativa importante. Esta es no contemplar el aprendizaje de las condiciones de normalidad  y anormalidades (potenciales ataques) del sistema más que en una etapa inicial de entrenamiento, imposibilitando la adaptación a nuevas condiciones normales variantes del sistema, a no ser que se realice un nuevo entrenamiento del sistema. Esto ocurre por el hecho que las técnicas tradicionales de aprendizaje automático no contemplan el problema del aprendizaje en tiempo real, impidiendo la continuidad operativa efectiva de estos sistemas en el mediano plazo sin detenciones por motivos de  actualización.

##Contexto Teórico

Se presenta una propuesta de sistema clasificador desarrollado en base a una población de agentes evolutivos. El modelo se centra en el efecto que producen cambios o perturbaciones del medio ambiente sobre individuos muy sensibles a estos, utilizándose el concepto de bioindicador [1] [2], esto es, la cuantificación de dicho efecto sobre la población de individuos para la determinación de anormalidades. Debido a que este medio ambiente es la representación continua de las características de un sistema monitoreado, el modelo presentado puede ser utilizado para la detección de anomalías sobre cualquier sistema caracterizable a partir de un flujo de parámetros que represente su estado.

La población de especies indicadoras artificiales (EIA) o bioindicadores artificiales logra su sensibilidad a cambios ambientales debido a un ajuste de los parámetros del modelo altamente desfavorable para la consolidación de la población de agentes. Con esto se logra un aumento de la presión de selección sobre los individuos de la población y una baja en su diversidad genética, lo cual según la evidencia experimental encontrada, amplifica los trastornos en la población ante perturbaciones, siendo de esta forma más fácil observar un cambio ante una anomalía.
La detección de anomalías es un enfoque de solución sobre el problema de clasificación, el cual consiste en la tarea de segregar objetos en un conjunto de clases diferentes. En algunos casos estas clases están predefinidas y no cambian en el tiempo, mientras en otros casos más complejos, las clases pueden no estar definidas a priori, y además estás pueden cambiar en el tiempo. Uno de estos escenarios complejos es el de Detección de Intrusos en Red (NIDS). En este dominio, el algoritmo clasificador debe enfrentar al menos dos clases fundamentales: tráfico normal y tráfico intrusivo. Estas clases no son estáticas y cambian por la variación habitual del comportamiento de los usuarios del sistema o por la presencia de un ataque nuevo/desconocido. Por esta razón se selecciona este escenario para probar las capacidades del clasificador propuesto en este trabajo.
El modelo EIA propone un enfoque ecológico, entendiendo que una población de agentes que se adapta plásticamente a los cambios del entorno para subsistir desarrolla tareas de aprendizaje, que en este contexto es su modificación estructural. Este enfoque ecológico se encuentra presente en la visión constructivista de Varela sobre el Sistema Inmune Biológico (BIS) [12], el cual pone énfasis en la autoafirmación y en el potencial homeostático, siendo esta la inspiración para Nanas [13] el cual implementa una red de términos adaptativa utilizada para el filtrado de información, que a diferencia del enfoque propuesto en este trabajo se basa en una red y nos en una población de agentes. 
La metáfora del sistema inmune ha sido ampliamente usada para la detección de intrusos en sistemas computacionales debido a que suponen objetivos similares: la detección y eliminación de agentes no propios/dañinos/desestabilizadores. Siendo justamente esta diferencia de conceptos la que ha dado pie a un prolífico y diverso [14] conjunto de técnicas híbridas nombrada en conjunto como Sistemas Inmunes Artificiales (AIS) [15] . Todas estas propuestas intentan rescatar las capacidades de identificación, eliminación de amenazas, tolerancia a fallos y adaptabilidad de los Sistemas Inmunes Biológicos (BIS), a través de una serie de propuestas tales como: Red Inmune Formal (FIN) [11] la cual se basa en procesos de apoptosis (programmed cell death) e immunización controlado por citokinas (Messenger proteins), Selección Clonal (CLONALG), que propone la proliferación de detectores hábiles para la detección de antígenos y la exploración de estos buscando mejorar la afinidad a través de hipermutación somática [17], Selección Negativa (LISYS)  [18] el cual basado en la maduración de linfocitos-T en el Thimus descarta aquellos detectores que producen auto-imunidad para producir tolerancia inmunológica [15] y modelos basados en la red inmune de Jerne (aiNet), los cuales exploran la afinidad de los detectores asociados como una red [20]. Existiendo evidencia que estas técnicas no abordan de forma consistente el cambio de la normalidad, además de basarse en modelos parciales y no plenamente aceptados [21]. Además Bersini también establece que las propuestas basadas en la tradicional concepción del sistema inmune como un ente defensivo implícita en las técnicas nombradas es incorrecta, obteniendo cada vez más oposición de biólogos [22]. Para Bersini el aporte real del modelo de BIS para la ingeniería recae en el concepto de endogenous doublé plasticity [23], el cual establece que el sistema se ajusta estructuralmente durante su funcionamiento de forma continua y plástica, integrando nuevos elementos como descartando antiguos, siendo este cambio controlado por su dinámica interna. Esto está basada en heurísticas simples, como lo son compensar elementos débiles, mantener la diversidad y suprimir la redundancia, esto es mantener equilibrios, a través de mecanismos ecológicos, lo cual es el sustento del presente trabajo.
El enfoque provisto por el modelo EIA, posee ventajas operativas importantes sobre el enfoque predominante para el desarrollo de Sistemas Detectores de Intrusos en Red (NIDS), el cual se centra principalmente en el uso de clasificadores basados en firmas o patrones preconocidos de ataque [25], que tienen las desventajas de depender de su constante actualización para resultar prácticos, además de ser inefectivos ante ataques desconocidos [25]. Esto último, es uno de los focos principales de interés en el desarrollo de AIS, los cuales son más cercanos a las propuestas basadas en la detección de anomalías, que si bien proveen una solución al problema de los ataques novedosos tienen el inconveniente de llevar asociados un incremento significativo en los falsos positivos [26].
El problema de clasificación de la detección de intrusos también ha sido abordado por diversas técnicas de Aprendizaje Supervisado (SL) tales como redes neuronales [30], máquinas de vectores de soporte [32], programación genética [29], lógica difusa [30] y Redes Bayesianas [31] entre otras. Para una revisión de desarrollos actuales se pueden consultar los siguientes trabajos actualizados [34,35]. El problema de las técnicas de SL, es que necesitan de procesos separados de entrenamiento en condiciones controladas, además de datos de entrenamiento rotulados y de calidad, logrando abordar en la parte dependencia de actualizaciones basadas en conocimiento experto; pero persistiendo en la incapacidad de detectar ataques desconocidos y en la falta de adaptación a cambios de escenario debido a su entrenamiento estático inicial. 
El presente trabajo corresponde a la continuación de  una prueba de concepto enmarcada en el desarrollo de un Sistema de Detección de Anomalías en Flujo (SDAF) [31], el cual está basado en dos hipótesis que son sustentadas por la evidencia encontrada durante su desarrollo, estas son:
•	hipótesis 1: Una ecología artificial de especies evolutivas genéticamente capaces, puede adaptarse a su medio ambiente descrito por un flujo de características, de tal forma de hacerse sensibles a sus trastornos de forma cuantificable. 
Esta hipótesis establece que es posible utilizar una técnica basada en nociones de vida artificial y sistemas ecológicos que logre determinar un SDAF primitivo, que sea capaz de detectar anomalías en flujos de datos. 
•	hipótesis 2: Un clasificador basado en nociones de sistemas ecológicos puede lograr ventajas operativas sobre otras técnicas de aprendizaje automático en el problema de clasificación en sistemas de detección de intrusos. Estos en términos de adaptación a cambios de las condiciones de normalidad y la detección de ataques innovadores. 

Ambas hipótesis fueron ciertas, siendo verificable las demostraciones en [35][36], por lo cual en este trabajo, corresponde llevar el concepto a un algoritmo computacionalmente eficiente para ser puesto en producción en sistemas de detección de intrusos reales.

##Referencias:
[1] Jeffrey, D.W., Madden, B.: Bioindicators and environmental management. Academic Press, London ; San Diego (1991).
[2] Gadzala-Kopciuch, R., Berecka, B., Bartoszewicz, J., Buszewski, B.: Some considerations about bioindicators in environmental monitoring. Polish Journal of Environmental Studies 13(5) (2004) 453{462 cited By (since
1996).
[3] Berrios, L., Prototipo de un sistema de prevención de intrusiones usando correlación de eventos de seguridad, in Departamento de Ingeniería Informática. 2007, Universidad de Santiago de Chile: Santiago, Chile.
[4] Varela, F.: El Fenómeno de la Vida. 2o edn. OCEANO, Santiago de Chile (2000) Nanas, N., de Roeck, A.: Autopoiesis, the immune system, and adaptive information ltering. Natural Computing 8 (2009) 387.
[5] Medina, C., Prototipo de Sistema de Detección de Intrusos usando Support Vector Machine, in Departamento de Ingeniería Informática. 2004, Universidad de Santiago de Chile: Santiago.
 [6] Coutinho,A.: A walk with francisco varela from rst- to second- generation networks: In search of the structure, dynamics and metadynamics of an organism-centered immune system. Biological Research 36(1) (2003) 17-26.
[7] Sistemas de detección de intrusos.  2004; Available from: http://escert.upc.edu/index.php/web/es/publicacion,178,3.html.
 [8] Estevez-Tapiador, J.M., P. Garcia-Teodoro, and J.E. Diaz-Verdejo, Anomaly detection methods in wired networks: a survey and taxonomy. Computer Communications, 2004. 27(16): p. 1569-1584.
[9] D. Stopel, R.M., Z. Boger, Y. Shahar, y Y. Elovici, Using artificial neural networks to detect unknown computer worms. Neural Computing & Applications, 2009. 18: p. 663-674.
 [10] Glickman, M., Balthrop, J., Forrest, S.: A machine learning evaluation of an artificial immune system. Evolutionary Computation 13(2) (2005) 179-212.
[11] Tarakanov, A.O.: Immunocomputing for intelligent intrusion detection. IEEE Computational Intelligence Magazine 3(2) (2008) 22-30.
[12] Cutello, V., Narzisi, G., Nicosia, G., Pavone, M.: Clonal selection algorithms: A comparative case study using e
ective mutation potentials. In Jacob, C and Pilat, ML and Bentley, PJ and Timmis, J, ed.: Articial Inmune Systems, Proceedings. Volume 3627 of Lecture Notes in Computer Science., Univ Calgary, Fac Sci; Univ Calgary, Dept Biochem; Univ Calgary, Dept Comp Sci; ARTIST; iCORE; MITACS; PIMS (2005) 13-28 4th International Conference on Artificial Immune Systems, Ban, Canada, Aug 14-17, 2005.
[13] Forrest, S., Perelson, A., Allen, L., Cherukuri, R.: Self-Nonself Discrimination in a Computer. In: 1994 IEEE Computer Society Symposium on Research in Security and Privacy, Proceedings, IEEE, COMP SOC; IEEE,
COMP SOC, TECH COMM SECUR & PRIVACY; INT ASSOC CRYPTOL RES (1994) 202{212 1994 IEEE-Computer-Society Symposium on Research in Security and Privacy, Oakland, CA, May 16-18, 1994.
[14] Harmer, P., Williams, P., Gunsch, G., Lamont, G.: An articial immune system architecture for computer security applications. Evolutionary Computation, IEEE Transactions on 6(3) (jun 2002) 252-280.
[15] de Castro, L., Von Zuben, F.: ainet an artificial immune network for data analysis. In Publishing, I.G., ed.: Data Mining: A Heuristic Approach. Idea Group Publishing (2001) 231-259.
[16] Wu, S.X., Banzhaf, W.: The use of computational intelligence in intrusion detection systems: A review. Applied Soft Computing 10(1) (2010) 1-35.
[17] Bersini, H.: Self-assertion versus self-recognition: A tribute to Francisco Varela. In Timmis, J., Bentley, P.J., eds.: Proceedings of the 1st International Conference on Artificial Immune Systems (ICARIS), University of Kent at Canterbury, University of Kent at Canterbury Printing Unit (September 2002) 107-112
[18] Dasgupta, D.: Artificial immune systems and their applications. Springer (1998).
[9] Mukherjee, B., Heberlein, L., Levitt, K.: Network intrusion detection. Network, IEEE 8(3) (MAY/JUN 1994) 26 -41.
[20] Greitzer, F.L., Moore, A.P., Cappelli, D.M., Andrews, D.H., Carroll, L.A., Hull, T.D.: Combating the insider cyber threat. IEEE Security & Privacy 6(1) (2008) 61-64.
[21] Estevez-Tapiador, J.M., Garcia-Teodoro, P., Diaz-Verdejo, J.E.: Anomaly detection methods in wired networks: a survey and taxonomy. Computer Communications 27(16) (2004) 1569-1584.
[22] Sklar, E.: Software review: NetLogo, a multi-agent simulation environment. Artificial Life 13(3) (SUM 2007) 303-311.
[23] Wilcoxon, F.: Individual Comparisons by Ranking Methods. Biometrics Bulletin 1(6) (1945) 80-83
[24] M., J., Halley: Ecology, evolution and 1f-noise. Trends in Ecology &amp;
Evolution 11(1) (1996) 33-37.
[25] Lippmann, R., Haines, J.W., Fried, D.J., Korba, J., Das, K.: The 1999 DARPA o
-line intrusion detection evaluation. Computer Networks the International Journal of Computer and Telecommunications Networking 34(4) (2000) 579-595
[26] Fawcett, T.: An introduction to ROC analysis. Pattern Recognition Letters 27(8) (2006) 861-874 ROC Analysis in Pattern Recognition.
[27] Olusola, A.A., Oladele, A.S., Abosede, D.O.: Analysis of KDD `99 Intrusion Detection Dataset for Selection of Relevance Features. In Ao, SI and Douglas, C and Grundfest, WS and Burgstone, J, ed.: World Congress om
Engineering and Computer Science, Vols 1 and 2. Lecture Notes in Engineering and Computer Science, Int Assoc Engn (2010) 162-168 World Congress on Engineering and Computer Science, San Francisco, CA, OCT
20-22, 2010.
[28] Tavallaee, M., Bagheri, E., Lu, W., Ghorbani, A.: A detailed analysis of the KDD CUP 99 data set. In: Computational Intelligence for Security and Defense Applications, 2009. CISDA 2009. IEEE Symposium on. (july 2009) 1-6.
[29] Kukielka, P., Kotulski, Z.: Analysis of Dierent Architectures of Neural Networks for Application in Intrusion Detection Systems. In Ganzha, M and Paprzycki, M and PelechPilichowski, T, ed.: 2008 International Mul-
ticonference on Computer Science and Information Technology (IMCSIT), Vols 1 and 2, IEEE (2008) 752{756 International Multiconference on Computer Science and Information Technology, Wisla, POLAND, OCT 20-22,
2008.
[30] Linda, O., Vollmer, T., Manic, M.: Neural Network Based Intrusion Detection System for Critical Infrastructures. In: IJCNN: 2009 International Joint Conference on Neural Networks, Vols 1- 6. IEEE International Joint Conference on Neural Networks (IJCNN), Int Neural Network Soc; IEEE Computat Intelligence Soc (2009) 102-109 International Joint Conference on Neural Networks, Atlanta, GA, JUN 14-19, 2009.
[31] S., H.: Neural Networks and Learning Machines. Third edition - new York edn. Prentice Hall (2009)
[32] Atreas, N., Karanikas, C., Tarakanov, A.: Signal processing by an immune type Tree Transform. In Timmis, J and Bentley, P and Hart, E, ed.: Artificial Inmune Systems, Proceedings. Volume 2787 of Lecture Notes in Com-
puter Science., EVONET; Hewlett Packard plc (2003) 111-119 2nd International Conference on Artificial Immune Sysyems, Edimburgh, Scotland, Sep 01-03, 2003.
[33] Horn, R., Johnson, C.: Matrix Analysis. Cambridge University Press (1986).
[34] Iturbe, J.A.: Modicación de los algoritmos basados en la metáfora del sistema inmunológico para la detección dinámica de intrusos en redes de datos. Master's thesis, Universidad de Santiago de Chile (2010).
[35] Pinacho Davidson, P.: Especies indicadoras artificiales (EIA), una técnicas de detección de anomalías en flujo basada en vida artificial. Master's thesis, Departamento de Ingeniería Informática, Universidad de Santiago de Chile (2011).
[36] Pedro Pinacho, Iván Pau, Max Chacón, Sergio Sánchez, An ecological approach to anomaly detection: the EIA model,ICARIS’12, Taormina, Italy, August 2012.

#SOBRE EL MODELO
Modelo EIA (AB) desarrollado por Pedro Pinacho Davidson
Estudiante de Doctorado de la EHU, bajo supervisión de doctores Lozano y Blum.


 
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="dos00-09" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>activeDelay2 != 0</exitCondition>
    <metric>activeDelay1</metric>
    <metric>activeDelay2</metric>
    <enumeratedValueSet variable="archivoAtaque">
      <value value="&quot;DoS-00&quot;"/>
      <value value="&quot;DoS-01&quot;"/>
      <value value="&quot;DoS-02&quot;"/>
      <value value="&quot;DoS-03&quot;"/>
      <value value="&quot;DoS-04&quot;"/>
      <value value="&quot;DoS-05&quot;"/>
      <value value="&quot;DoS-06&quot;"/>
      <value value="&quot;DoS-07&quot;"/>
      <value value="&quot;DoS-08&quot;"/>
      <value value="&quot;DoS-09&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="activar-penalizacion">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Gold00-09" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>activeDelay2 != 0</exitCondition>
    <metric>ticks - 950activeDelay1</metric>
    <metric>activeDelay2</metric>
    <enumeratedValueSet variable="archivoAtaque">
      <value value="&quot;Gold-00&quot;"/>
      <value value="&quot;Gold-01&quot;"/>
      <value value="&quot;Gold-02&quot;"/>
      <value value="&quot;Gold-03&quot;"/>
      <value value="&quot;Gold-04&quot;"/>
      <value value="&quot;Gold-05&quot;"/>
      <value value="&quot;Gold-06&quot;"/>
      <value value="&quot;Gold-07&quot;"/>
      <value value="&quot;Gold-08&quot;"/>
      <value value="&quot;Gold-09&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="activar-penalizacion">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="LOIC00-09" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>activeDelay2 != 0</exitCondition>
    <metric>activeDelay1</metric>
    <metric>activeDelay2</metric>
    <enumeratedValueSet variable="archivoAtaque">
      <value value="&quot;LOIC-00&quot;"/>
      <value value="&quot;LOIC-01&quot;"/>
      <value value="&quot;LOIC-02&quot;"/>
      <value value="&quot;LOIC-03&quot;"/>
      <value value="&quot;LOIC-04&quot;"/>
      <value value="&quot;LOIC-05&quot;"/>
      <value value="&quot;LOIC-06&quot;"/>
      <value value="&quot;LOIC-07&quot;"/>
      <value value="&quot;LOIC-08&quot;"/>
      <value value="&quot;LOIC-09&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="activar-penalizacion">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Tor00-09" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>activeDelay2 != 0</exitCondition>
    <metric>activeDelay1</metric>
    <metric>activeDelay2</metric>
    <enumeratedValueSet variable="archivoAtaque">
      <value value="&quot;Tor-00&quot;"/>
      <value value="&quot;Tor-01&quot;"/>
      <value value="&quot;Tor-02&quot;"/>
      <value value="&quot;Tor-03&quot;"/>
      <value value="&quot;Tor-04&quot;"/>
      <value value="&quot;Tor-05&quot;"/>
      <value value="&quot;Tor-06&quot;"/>
      <value value="&quot;Tor-07&quot;"/>
      <value value="&quot;Tor-08&quot;"/>
      <value value="&quot;Tor-09&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="activar-penalizacion">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Switch00-09" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>activeDelay2 != 0</exitCondition>
    <metric>activeDelay1</metric>
    <metric>activeDelay2</metric>
    <enumeratedValueSet variable="archivoAtaque">
      <value value="&quot;Switch-00&quot;"/>
      <value value="&quot;Switch-01&quot;"/>
      <value value="&quot;Switch-02&quot;"/>
      <value value="&quot;Switch-03&quot;"/>
      <value value="&quot;Switch-04&quot;"/>
      <value value="&quot;Switch-05&quot;"/>
      <value value="&quot;Switch-06&quot;"/>
      <value value="&quot;Switch-07&quot;"/>
      <value value="&quot;Switch-08&quot;"/>
      <value value="&quot;Switch-09&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="activar-penalizacion">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="LOIC" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="archivoAtaque">
      <value value="&quot;LOIC-00&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="activar-penalizacion">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
