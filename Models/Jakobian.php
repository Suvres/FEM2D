<?php

declare(strict_types=1);

namespace Mes\Models;

use JetBrains\PhpStorm\Pure;
use Mes\Data\Data;

class Jakobian {

    /**
     * @var float[][]
     */
    public array $jakobian = [
        [0.0, 0.0],
        [0.0, 0.0]
    ];

    /**
     * @var float[][]
     */
    public array $jakobian_inv = [
        [0.0, 0.0],
        [0.0, 0.0]
    ];

    /**
     * @var float
     */
    public float $det_jakobian = 0.0;

    // generowanie jakobianów dla elementów i ich punktów całkowania
    public static function generateJakobian(int $id, Element4_2D $element4_2D, Grid $grid): void
    {
        for($i = 0; $i < Data::DATA_points**2; $i++) {
            /** @var Element $element */
            $element = $grid->elements[$id];

            $element->jakobiany[] = self::jakobian($i, $id, $element4_2D, $grid);
        }
    }

    private static function jakobian(int $pc, int $id, Element4_2D $element4_2D, Grid $grid): Jakobian
    {
        /** @var Element $element */
        $element = $grid->elements[$id];

        $jac = new self();
        //dX / dKsi: czyli suma (dN_i/dKsi) * x_i
        $jac->jakobian[0][0] =
            $element4_2D->dN_dKsi[$pc][0] * $grid->nodes[$element->ID[0]]->x
            + $element4_2D->dN_dKsi[$pc][1] * $grid->nodes[$element->ID[1]]->x
            + $element4_2D->dN_dKsi[$pc][2] * $grid->nodes[$element->ID[2]]->x
            + $element4_2D->dN_dKsi[$pc][3] * $grid->nodes[$element->ID[3]]->x;

        //dX / dEta: czyli suma (dN_i/dEta) * x_i
        $jac->jakobian[0][1] =
            $element4_2D->dN_dEta[$pc][0] * $grid->nodes[$element->ID[0]]->x
            + $element4_2D->dN_dEta[$pc][1] * $grid->nodes[$element->ID[1]]->x
            + $element4_2D->dN_dEta[$pc][2] * $grid->nodes[$element->ID[2]]->x
            + $element4_2D->dN_dEta[$pc][3] * $grid->nodes[$element->ID[3]]->x;

        //dY / dKsi: czyli suma (dN_i/dKsi) * y_i
        $jac->jakobian[1][0] =
            $element4_2D->dN_dKsi[$pc][0] * $grid->nodes[$element->ID[0]]->y
            + $element4_2D->dN_dKsi[$pc][1] * $grid->nodes[$element->ID[1]]->y
            + $element4_2D->dN_dKsi[$pc][2] * $grid->nodes[$element->ID[2]]->y
            + $element4_2D->dN_dKsi[$pc][3] * $grid->nodes[$element->ID[3]]->y;

        //dY / dEta: czyli suma (dN_i/dEta) * y_i
        $jac->jakobian[1][1] =
            $element4_2D->dN_dEta[$pc][0] * $grid->nodes[$element->ID[0]]->y
            + $element4_2D->dN_dEta[$pc][1] * $grid->nodes[$element->ID[1]]->y
            + $element4_2D->dN_dEta[$pc][2] * $grid->nodes[$element->ID[2]]->y
            + $element4_2D->dN_dEta[$pc][3] * $grid->nodes[$element->ID[3]]->y;

        //wyznacznik ze wzoru na wyznacznik macierzy
        $jac->det_jakobian = ($jac->jakobian[0][0] * $jac->jakobian[1][1]) - ($jac->jakobian[1][0] * $jac->jakobian[0][1]);

        // [jakobian] do ^-1
        $jac->jakobian_inv[0][0] = $jac->jakobian[1][1] / $jac->det_jakobian;
        $jac->jakobian_inv[1][1] = $jac->jakobian[0][0] / $jac->det_jakobian;
        $jac->jakobian_inv[1][0] = - $jac->jakobian[1][0] / $jac->det_jakobian;
        $jac->jakobian_inv[0][1] = - $jac->jakobian[0][1] / $jac->det_jakobian;

        return $jac;
    }

}
