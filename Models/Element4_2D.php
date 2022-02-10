<?php

declare(strict_types=1);

namespace Mes\Models;

use Mes\Data\Data;

class Element4_2D
{
    public array $dN_dKsi = [];
    public array $dN_dEta = [];
    public array $gauss = [];

    public function __construct() {
        $this->gauss = (new GaussTable(Data::DATA_points))->quadrat;
        $c = 0;
        for($i = 0; $i < Data::DATA_points; $i++) {
            for($j = 0; $j < Data::DATA_points; $j++) {
                $this->dN_dKsi[$c][0] = -(1.0 / 4.0) * (1.0 - $this->gauss[0][$i]);
                $this->dN_dKsi[$c][1] = (1.0 / 4.0) * (1.0 - $this->gauss[0][$i]);
                $this->dN_dKsi[$c][2] = (1.0 / 4.0) * (1.0 + $this->gauss[0][$i]);
                $this->dN_dKsi[$c][3] = -(1.0 / 4.0) * (1.0 + $this->gauss[0][$i]);

                $this->dN_dEta[$c][0] = -(1.0 / 4.0) * (1.0 - $this->gauss[0][$j]);
                $this->dN_dEta[$c][1] = -(1.0 / 4.0) * (1.0 + $this->gauss[0][$j]);
                $this->dN_dEta[$c][2] = (1.0 / 4.0) * (1.0 + $this->gauss[0][$j]);
                $this->dN_dEta[$c][3] = (1.0 / 4.0) * (1.0 - $this->gauss[0][$j]);
                $c++;
            }
        }

    }

}