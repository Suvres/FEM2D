<?php

declare(strict_types=1);

namespace Mes\Models;

class Node{
    public function __construct(
        public float $x,
        public float $y,
        public float $t0,
        public float $bc = 0
    ){}
}

