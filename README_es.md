# Herramienta de ingeniería y modelado de radiofrecuencias

<!-- hy-mt2-i18n:start -->
[English](./README.md) | [中文](./README_zh-CN.md) | [日本語](./README_ja.md) | **Español**
<!-- hy-mt2-i18n:end -->


Este proyecto ofrece un marco de simulación de alta fidelidad para sistemas de radar, telecomunicaciones y otros sistemas de ingeniería de radiofrecuencias. Su objetivo es lograr una modelización y simulación realistas a nivel de diseño, siguiendo los estándares y metodologías presentes en la literatura técnica e ingenieril avanzada.

## Características principales

- **Modelado integral del sistema:**  
  Modela los componentes clave de los sistemas de radar y telecomunicaciones, incluyendo cálculos analíticos de la sección transversal de radar (RCS) y patrones de radiación de antenas.  
- **Núcleos optimizados para el rendimiento:**  
  Implementa algoritmos altamente optimizados que aprovechan las instrucciones intrínsecas de Intel (SSE/AVX/AVX2/AVX512) para una vectorización manual masiva. Se utiliza la autovectorización a nivel de compilador para estadísticas descriptivas y análisis de rendimiento.  
- **Aceleración con GPGPU:**  
  Incluye una amplia base de código CUDA (~15,000 líneas) que abarca núcleos computacionales y rutinas auxiliares para simulaciones aceleradas por GPU.  
- **Arquitectura modular:**  
  Está organizada como un conjunto de módulos independientes, cada uno de los cuales describe componentes modelados específicos. La biblioteca puede funcionar como backend computacional o integrarse con una interfaz gráfica de usuario.  
- **Alcance de los componentes:**  
  El marco está estructurado en torno a cuatro dominios principales de simulación:  
  1. Modelado y simulación de sistemas de telecomunicaciones en RF (parcialmente implementado).  
  2. Modelado y simulación de radares (parcialmente implementado).  
  3. Procesamiento de señales en banda base (envoltura compleja) y modelado de receptores óptimos (aún no implementado).  
  4. Bibliotecas de soporte (por ejemplo, atmósfera, terreno y edificios) (parcialmente implementado).
  
## Visión general de la implementación

- **Rutas de ejecución SIMD:**  
  Cientos de núcleos computacionales implementados para precisión doble y simple, centrados en el modelado analítico de RCS y antenas.  
- **Ruta CUDA:**  
  Núcleos para GPU que permiten cálculos a gran escala y de alto rendimiento.  
- **Estado actual:**  
  El proyecto cuenta con cientos de núcleos SIMD y un conjunto completo de herramientas analíticas y de simulación para el análisis y modelado de sistemas de radar y telecomunicaciones.

## Uso

Este software está diseñado como una biblioteca de cálculo de backend para aplicaciones avanzadas de simulación y modelado. Puede integrarse en entornos de software más amplios o conectarse a interfaces gráficas de usuario con el fin de realizar visualizaciones y análisis.

## Contribuciones

Se dan la bienvenida las contribuciones, especialmente de quienes cuenten con experiencia en:  
- Sistemas de radar, telecomunicaciones e ingeniería de radio  
- Computación de alto rendimiento (SIMD, CUDA)  
- Métodos numéricos y computación científica

Si está interesado, por favor cree un problema o una solicitud de pull request.

## Licencia

Este proyecto está licenciado bajo los términos de la licencia GPLv3.

## Agradecimientos

Este proyecto se basa en la literatura de ingeniería y técnica, y está impulsado por el compromiso con un modelado de sistemas realista y eficiente.






