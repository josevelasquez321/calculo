# Simulador del sistema de frenado · Proyecto ABPro

Página web interactiva del proyecto **«Optimización cinemática y dinámica del sistema de frenado mediante modelado matemático y cálculo diferencial»** (Cálculo Diferencial, INACAP), aplicado a un Chevrolet Sail.

**Ver en línea:** https://josevelasquez321.github.io/calculo/

## Pestañas

1. **Velocidad y frenado**: x(t), v(t) = dx/dt, a(t) = dv/dt, distancia de frenado y de detención, animación de la frenada y comparación entre dos velocidades.
2. **Masa y carga**: efecto de la masa (límite de los frenos frente al límite del neumático) y transferencia de carga entre ejes.
3. **Energía y temperatura**: energía cinética, dT/dt durante la frenada y temperatura del disco en frenadas repetidas.
4. **Presión hidráulica**: P(t), dP/dt y distancia extra por la subida de presión (manguera de goma frente a línea acerada).

La página es un solo archivo (`index.html`) sin dependencias externas. La versión en Python que genera los gráficos del informe está en `python/simulador_frenado.py`.

Las pruebas del proyecto se realizan en esta página, no en terreno. Los parámetros marcados como supuestos (masa del disco, constantes de tiempo, altura del centro de gravedad) son valores de referencia.
