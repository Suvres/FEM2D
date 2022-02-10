<?php

declare(strict_types=1);

namespace Mes\Models;

class GaussTable {
    public array $quadrat = [];

    public function __construct(int $i)
    {
        $quadrat = [
            2 => [
                [ -1.0 / sqrt(3), 1.0 / sqrt(3) ],
                [1,1]
            ],
            3 => [
                [-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)],
                [5.0/9.0, 8.0/9.0, 5.0/9.0]
            ],
            4 => [
                [
                    -sqrt((3.0/7.0) + ((2.0/7.0)*6.0/5.0)),
                    -sqrt((3.0/7.0) - ((2.0/7.0)*6.0/5.0)),
                    sqrt((3.0/7.0) - ((2.0/7.0)*6.0/5.0)),
                    sqrt((3.0/7.0) + ((2.0/7.0)*6.0/5.0))],
                [
                    (18.0 - sqrt(30.0) / 36.0),
                    (18.0 + sqrt(30.0) / 36.0),
                    (18.0 + sqrt(30.0) / 36.0),
                    (18.0 - sqrt(30.0) / 36.0)
                ]
            ],
        ];

        $this->quadrat = $quadrat[$i];
    }
}
