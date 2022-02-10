<?php

declare(strict_types=1);

namespace Mes\Models;

use Mes\Data\Data;
use Mes\Service\Helper;

class Element {
    /**
     * @var float[][]
     */
    public array $h;

    /**
     * @var float[][]
     */
    public array $hbc;

    /**
     * @var float[]
     */
    public array $p;

    /**
     * @var float[][]
     */
    public array $c;

    /**
     * @var Jakobian[]
     */
    public array $jakobiany;

    /**
     * @var float[][]
     */
    private array $dn_x;

    /**
     * @var float[][]
     */
    private array $dn_y;

    /**
     * @var array{ID1: int, ID2: int, ID3: int, ID4: int}
     */
    public function __construct(
        /**
         * @var float[][]
         */
        public array $ID,
    ){
        for ($i = 0; $i < 4; $i++) {
            for($j = 0; $j < 4; $j++) {
                $this->hbc[$i][$j] = 0;
                $this->c[$i][$j] = 0;
            }

            $this->p[$i] = 0;
        }
        for($i = 0; $i < 4; $i++) {
            $this->h[$i] = [0,0,0,0];
        }
    }

    public function  computeHmatrix(Element4_2D $element4_2D): void
    {
        for($j = 0; $j < Data::DATA_points**2; $j++) {
            [$jac_inv_0, $jac_inv_1] = $this->jakobiany[$j]->jakobian_inv;

            // Wyliczenie dN / dX dla każdej funkcji kształtu
            $dn_dX = [
                ( $jac_inv_0[0] * $element4_2D->dN_dKsi[$j][0] ) + ( $jac_inv_0[1] * $element4_2D->dN_dEta[$j][0] ),
                ( $jac_inv_0[0] * $element4_2D->dN_dKsi[$j][1] ) + ( $jac_inv_0[1] * $element4_2D->dN_dEta[$j][1] ),
                ( $jac_inv_0[0] * $element4_2D->dN_dKsi[$j][2] ) + ( $jac_inv_0[1] * $element4_2D->dN_dEta[$j][2] ),
                ( $jac_inv_0[0] * $element4_2D->dN_dKsi[$j][3] ) + ( $jac_inv_0[1] * $element4_2D->dN_dEta[$j][3] )
            ];

            // Wyliczenie dN / dY dla każdej funkcji kształtu
            $dn_dY = [
                ( $jac_inv_1[0] * $element4_2D->dN_dKsi[$j][0] ) + ( $jac_inv_1[1] * $element4_2D->dN_dEta[$j][0] ),
                ( $jac_inv_1[0] * $element4_2D->dN_dKsi[$j][1] ) + ( $jac_inv_1[1] * $element4_2D->dN_dEta[$j][1] ),
                ( $jac_inv_1[0] * $element4_2D->dN_dKsi[$j][2] ) + ( $jac_inv_1[1] * $element4_2D->dN_dEta[$j][2] ),
                ( $jac_inv_1[0] * $element4_2D->dN_dKsi[$j][3] ) + ( $jac_inv_1[1] * $element4_2D->dN_dEta[$j][3] )
            ];

//                Helper::printmatrix($dn_dX);
//                Helper::printmatrix($dn_dY);
//                printf("--- %s ---\n", $j);

            /**
             * modulo — po kolumnach przykład: 4 mod 3 = 1
             * dzielenie — po wierszach przykład: przedział [4, 6] / 3 = 1
             */
            $this->generateHMatrix2(
                ($element4_2D->gauss[1][$j % Data::DATA_points] *
                    $element4_2D->gauss[1][(int)($j / Data::DATA_points)]),
                $j, $dn_dX, $dn_dY);
        }
    }

    private function generateHMatrix2(float $gausWeight ,int $point, array $dn_X, array $dn_Y): void
    {
        for($j = 0; $j < 4; $j++) {
            for($k = 0; $k < 4; $k++) {
                $this->h[$j][$k] = $this->h[$j][$k] ?? 0.0;

                $this->h[$j][$k] +=
                    (($dn_X[$j] * $dn_X[$k]) + ($dn_Y[$j] * $dn_Y[$k]))
                    * $gausWeight * Data::STALA_K
                    * $this->jakobiany[$point]->det_jakobian;
            }
        }
    }

    // Wyliczenie warunku brzegowego dla elementu skończonego
    // czyli macierz Hbc i macierz P
    public function computeHbcAndPMatrix(Element4_2D $element4_2D, Grid $grid): void
    {
        // przechodzi po nodach
        for($i = 0; $i < 4; $i++) {
            $i_2 = $i + 1 > 3 ? 0 : $i + 1;

            $node = $grid->nodes[$this->ID[$i]];
            $node1 = $grid->nodes[$this->ID[$i_2]];

            \assert($node instanceof Node);
            \assert($node1 instanceof Node);

            // sprawdzenie, czy ściana może mieć obliczony warunek brzegowy,
            // czyli czy ściana jest krańcowa
            if($node->bc === 0.0 || $node1->bc === 0.0) {
                continue;
            }

            // Długość odcinka w ukłądzie kartezjańskiom
            $det_jac = sqrt((($node->x - $node1->x) ** 2) + (($node->y - $node1->y) ** 2)) / 2.0;

            // wyliczenie Hbc i P dla każdej ściany -> ten for odpowiada za ścianę, bo
            // w niej są punkty całkowania, całkowanie po powierzchni
            for($j = 0; $j < Data::DATA_points; $j++) {
                $compute_functions = [0.0, 0.0, 0.0, 0.0];

                //Znak funkcji kształtu, trochę pokombinowane, ale ułatwia korzystanie ze ścian
                $Ksi_znak = [-1.0, 1.0, 1.0, -1.0];
                $Eta_znak = [-1.0, -1.0, 1.0, 1.0];

                // Funkcja Kształtu N = 1/4 * (1 - Ksi) * (1 - Eta)
                $Ksi = $i === 1 || $i === 3 ? 1 : $element4_2D->gauss[0][$j];
                $Eta = $i === 0 || $i === 2 ? 1 : $element4_2D->gauss[0][$j];

                // Ten if mimo, że wygląda dziwnie, bo obsługuje schemat całkowania 3pkt,
                // czyli pozwala na to by ksi albo eta była równa 0, na środku ściany
                if( Data::DATA_points % 3 === 0) {
                    $Ksi = ($i === 0 || $i === 2 ) && $j === (int) floor(Data::DATA_points / 2) ? 0.0 : $Ksi;
                    $Eta = ($i === 1 || $i === 3 ) && $j === (int) floor(Data::DATA_points / 2) ? 0.0 : $Eta;
                }

                // Odwracanie znaków na ścianach: rozwiązanie problemu siatki kartezjańskiej,
                // czyli po prostu jak mamy siatkę, to w zależności od pozycji na układzie współrzędnych
                // mamy inne znaki
                $Ksi = $i >= 2 ? $Ksi * -1 : $Ksi;
                $Eta = $i === 3 || $i === 0 ? $Eta * -1 : $Eta;

                // wyliczanie funkcji kształtów potrzebnych do obliczenia macierzy Hbc i P,
                // liczymy dwie funkcje kształtu w zależności od nodów, które są na ścianie
                $compute_functions[$i] = (1.0/4.0)
                    * (1.0 + ($Ksi_znak[$i] * $Ksi))
                    * (1.0 + ($Eta_znak[$i] * $Eta));

                $compute_functions[$i_2] = (1.0/4.0)
                    * (1.0 + ($Ksi_znak[$i_2] * $Ksi))
                    * (1.0 + ($Eta_znak[$i_2] * $Eta));

                // Liczenie macierzy Hbc, tutaj jest wzór integral_S = alfa({N} * {N}_trans)*dS
                // dS to det_j to długość odcinka / 2
                $this->computeHbc($element4_2D->gauss[1][$j], $det_jac, $compute_functions);

                // liczenie wektora P: integral_S =  alfa * {N} * t_0 * dS
                // tak jak poprzednio dS to det_j czyli długość odcinka / w
                $this->computePmatrix($j, $element4_2D->gauss[1][$j], $compute_functions, $det_jac);
            }
        }
    }

    // Macierz C - czyli pojemność cieplna elementu, liczona jest tak jak macierz H dla funkcji kształtu
    // każdego elementu ze wzoru integral_V = RO * C * ({N}*{N}_trans)dV
    // przy wyliczeniu tego det_j jest to wyznacznik z obiektu jakobian, który jest liczony wyżej
    public function computeCMatrix(Element4_2D $element4_2D): void
    {
        $c = 0;
        for($j = 0; $j < Data::DATA_points; $j++) {
            $pc1 = $element4_2D->gauss[0][$j];

            for($k = 0; $k < Data::DATA_points; $k++) {
                $pc2 = $element4_2D->gauss[0][$k];

                $compute_shapes = [
                    (1.0/4.0) * (1.0 - ($pc1)) * (1.0 - ( $pc2)),
                    (1.0/4.0) * (1.0 + ($pc1)) * (1.0 - ( $pc2)),
                    (1.0/4.0) * (1.0 + ($pc1)) * (1.0 + ( $pc2)),
                    (1.0/4.0) * (1.0 - ($pc1)) * (1.0 + ( $pc2))
                ];

                $this->generateCMatrix(
                    $element4_2D->gauss[1][$j] * $element4_2D->gauss[1][$k],
                    $c, $compute_shapes);
                $c++;
            }
        }
    }

    private function generateCMatrix(float $gausWeight ,int $point, array $compute_shapes): void
    {
        for($j = 0; $j < 4; $j++) {
            for($k = 0; $k < 4; $k++) {
                $this->c[$j][$k] +=
                    $compute_shapes[$j] * $compute_shapes[$k]
                    * Data::RO * Data::C * $gausWeight
                    * $this->jakobiany[$point]->det_jakobian;
            }
        }
    }

    /**
     * @param float[] $compute_shapes
     */
    private function computePmatrix(int $wall, float $weight, array $compute_shapes, float $det_jac): void
    {
        for($i = 0; $i < 4; $i++ ) {
            $this->p[$i] +=  $compute_shapes[$i] * Data::ALFA * $det_jac * Data::STALA_T_1[$wall] * $weight;
        }
    }

    private function computeHbc(float $weight, float $det_jac, array $shape_functions): void
    {
        $c = count($shape_functions);
        for($j = 0; $j < $c; $j++) {
            $this->hbc[0][$j] += $shape_functions[0] * $shape_functions[$j] * Data::ALFA * $weight * $det_jac;
            $this->hbc[1][$j] += $shape_functions[1] * $shape_functions[$j] * Data::ALFA * $weight * $det_jac;
            $this->hbc[2][$j] += $shape_functions[2] * $shape_functions[$j] * Data::ALFA * $weight * $det_jac;
            $this->hbc[3][$j] += $shape_functions[3] * $shape_functions[$j] * Data::ALFA * $weight * $det_jac;
        }
    }
}
