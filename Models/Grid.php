<?php

declare(strict_types=1);

namespace Mes\Models;

use Mes\Data\Data;

class Grid
{
    public float $H;
    public float $B;
    public float $nH;
    public float $nB;
    public int $nN;
    public int $nE;
    public array $nodes;
    public array $elements;

    public function __construct(float $H, float $B, float $nH, float $nB)
    {
        $this->H = $H;
        $this->B = $B;
        $this->nH = $nH;
        $this->nB = $nB;
        $this->nN = (int)($this->nH * $this->nB);
        $this->nE = (int)(($this->nH - 1) * ($this->nB - 1));
    }


    public static function generateGrid(float $H, float $B, float $nH, float $nB): grid
    {
        $grid = new Grid($H, $B, $nH, $nB);
        $nodes = [];
        $elements = [];

        $x = 0;
        $dY = ($grid->H / ($grid->nH - 1.0));
        $dX = ($grid->B / ($grid->nB - 1.0));

        for ($i = 0; $i < ($grid->nN); $i++ ) {
            $x = $i  < ($grid->nH + ($x * $grid->nH)) ? $x : $x+1;

            $x_tmp = $x * $dX;
            $y = ($i % $grid->nH) * $dY;

            $bc = $x_tmp === 0.0 || $y === 0.0 || $x_tmp === $grid->B || $y === $grid->H ? 1 : 0;

            $nodes[$i] = new Node(
                $x_tmp,
                $y,
                Data::STALA_T_0,
                $bc
            );
        }

        $e1 = 0;
        for ($i = 0; $i < $grid->nE; $i++) {
            $e1 = $i > 0 && $i % ($grid->nH - 1.0)  === 0 ? $e1 + 1 : $e1;

            $elements[$i] = new Element([
                $i + $e1,
                $i + $e1 + $grid->nH,
                $i + $e1 + $grid->nH + 1,
                $i + $e1 + 1,
            ]);
        }

        $grid->nodes = $nodes;
        $grid->elements = $elements;

        unset($elements, $nodes);
        return $grid;
    }

}
