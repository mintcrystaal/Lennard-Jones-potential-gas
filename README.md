# bMT_mc4statemodel

## Описание

Программа для моделирования частиц одноатомного и двухатомного газа с взаимодействием по потенциалу Леннарда-Джонса.

## Структура проекта

* `gas.py` - основной файл программы, моделирование одноатомного газа по принципу молекулярного моделирования
* 'diatomic_gas.py' - моделирование двухатомного газа по принципу молекулярного моделирования
* 'make_graphs.py' - программа, строящая графики для выполнения закона сохранения энергии
* 'brownian.py' - программа, строящая график для проверки выполнения закона Эйнштейна-Смолуховского


С помощью программы был вычислен коэффициент самодиффузии для различных одноатомных газов.


| Вещество | $\sigma$, Å | $\varepsilon/k$, К | $m \cdot 10^{-27}$ кг | $D$, $\frac{\text{см}^2}{\text{с}}$ | $D_{\text{table}}$, $\frac{\text{см}^2}{\text{с}}$ |
|----------|----------------|--------------------|------------------------|-----------------------------|----------------------------------|
| He       | 2.628          | 5.465              | 6.4                    | 0.80 $\pm$ 0.06             | 0.36                             |
| Ar       | 3.401          | 116.81             | 63.92                  | 1.51 $\pm$ 0.05             | 1.58                             |
| Ne       | 2.775          | 36.831             | 32.29                  | 0.97 $\pm$ 0.07             | 0.45                             |


