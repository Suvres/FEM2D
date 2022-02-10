<?php

declare(strict_types=1);

namespace Mes\Models;

use Mes\Data\Data;
use Mes\Service\Helper;

class MatrixHg {
    /**
     * @var float[][]
     */
    public static array $H;

    public static array $C;


    /**
     * @var float[]
     */
    public static array $P;

    /**
     * @var float[][]
     */
    public static array $matrixExpand;

    public static function init(int $nodes): void
    {
        for ($i = 0; $i < $nodes; $i++) {
            for($j = 0; $j < $nodes; $j++) {
                self::$H[$i][$j] = 0;
                self::$C[$i][$j] = 0;
            }

            self::$P[$i] = 0;
        }
    }

    // w skrócie, macierz rozszczerzona to [H] z doklejonym {P} czyli:
    // |- H_11 H_12 H_13 -|    { P_1 }    |- H_11 H_12 H_13 P_1-|
    // |  H_21 H_22 H_23  |  i { P_2 } => |  H_21 H_22 H_23 P_2 |
    // |- H_31 H_32 H_33 -|    { P_3 }    |- H_31 H_32 H_33 P_3-|
    public static function createExpandMatrix(): void
    {
        self::$matrixExpand = self::$H;
        $last = count(self::$H);

        $_counter = 0;
        foreach (self::$P as $p) {
            self::$matrixExpand[$_counter][$last] = $p;
            $_counter++;
        }
    }

    public static function generateMatrixHAndP(Grid $grid): void
    {
        $element4_2D = new Element4_2D();
        for($i = 0; $i < $grid->nE; $i++ ) {

            $element = $grid->elements[$i];
            \assert($element instanceof Element);

            // generowanie tablicy jakobianów
            // każdy element ma tablice jakobianów o wielkości Points**2
            // każdy punkt całkowania dla elementu 2D ma swój jakobian
            // wiążę się to z tym, żę jakobian to stosunek pól powierzchni układu globalnego i lokalnego
            Jakobian::generateJakobian($i, $element4_2D, $grid);

            // generowanie macierzy H dla każdego elementu wewnątrz niego
            $element->computeHmatrix($element4_2D);

            // generowanie macierzy Hbc i P czyli dwóch części warunku brzegowego, rozdzielonych ze względu na to,
            // że wzór na konwekcje to: q = a(t - t_0) tak więc mamy q = at - at_0
            // z tego znamy wszystko a, t_0, dlatego jedna część idzie do Hbc, a ta znana do P
            $element->computeHbcAndPMatrix($element4_2D, $grid);

            // macierz C służy do wyliczenia pojemności cieplnej materiału, dzięki czemu można obliczyć
            // temperaturę po kroku czasowym delta_T
            $element->computeCMatrix($element4_2D);

            // Agregacja wykorzystuje wzory:
            // [H] = [H] + [C]/delta_T
            // {P} + ({[C]/delta_T}*{T_0})
            self::aggregateHAndPMatrix($element, $grid);
        }

    }

    private static function aggregateHAndPMatrix(Element $element, Grid $grid): void
    {

        // agregacja macierzy H wykorzystując lokalne macierze H z każdego elementu
        // macierz Hbc (kawałek warunku brzegowego) oraz macierz C czyli pojemność cieplną
        for($j = 0; $j < 4; $j++) {
            for($k = 0; $k < 4; $k++) {
                self::$H[$element->ID[$j]][$element->ID[$k]] +=
                    $element->h[$j][$k]
                    + $element->hbc[$j][$k]
                    + ( $element->c[$j][$k] / Data::DELTA_T );

                // wersja do testów niezwiązana z rozwiązaniem, tylko by sprawdzić poprawność macierzy C
                self::$C[$element->ID[$j]][$element->ID[$k]] += $element->c[$j][$k];
            }
        }

        // agregacja wektora P wykorzystując wektory P lokalne dla elementu oraz macierz C
        for($j = 0; $j < 4; $j++) {
            $c = 0.0;

            // suma "wiersza" macierzy C dla kroku czasowego razy t_0
            for($k = 0; $k < 4; $k++) {
                $c += $element->c[$j][$k] / Data::DELTA_T
                    * $grid->nodes[$element->ID[$k]]->t0;
            }

            self::$P[$element->ID[$j]] += $element->p[$j] + $c;
        }
    }
}