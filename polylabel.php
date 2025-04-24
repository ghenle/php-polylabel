<?php
/**
 * PHP port of https://github.com/mapbox/polylabel
 * for use with https://gasparesganga.com/labs/php-shapefile/
 **/
namespace utility;

# prevent direct access to included files
(__FILE__ != $_SERVER['SCRIPT_FILENAME']) or exit('No');


class polylabel {
  CONST XKEY = 0;
  CONST YKEY = 1;

  // pole of inaccessibility
  public $centerOfMass;
  public $centroid;
  public $boundsCentroid;
  public $bounds;
  public $area;

  public function __construct($polygon, $precision = 1.0, $debugCallback = null) {

    // find the bounding box of the outer ring
    $this->bounds = array(
      'xmax' => \config::PHP_INT_MIN,
      'xmin' => \config::PHP_INT_MAX,
      'ymax' => \config::PHP_INT_MIN,
      'ymin' => \config::PHP_INT_MAX
    );

    foreach ( $polygon as $ring ) {
      foreach ( $ring as $xy ) {
        $this->bounds['xmin'] = min( $this->bounds['xmin'], $xy['x'] );
        $this->bounds['xmax'] = max( $this->bounds['xmax'], $xy['x'] );
        $this->bounds['ymin'] = min( $this->bounds['ymin'], $xy['y'] );
        $this->bounds['ymax'] = max( $this->bounds['ymax'], $xy['y'] );
      }
    }

    $width    = $this->bounds['xmax'] - $this->bounds['xmin'];
    $height   = $this->bounds['ymax'] - $this->bounds['ymin'];
    $cellSize = min($width, $height);

    if( $cellSize === 0 ) {
      return [
        'x' => $this->bounds['xmin'],
        'y' => $this->bounds['ymin'],
        'd' => 0
      ];
    }

    // a priority queue of cells in order of their "potential" (max distance to polygon)
    $cellQueue = new \utility\CellQueue();

    // cover polygon with initial cells
    $h = $cellSize / 2;
    for( $x = $this->bounds['xmin']; $x < $this->bounds['xmax']; $x += $cellSize ) {
      for( $y = $this->bounds['ymin']; $y < $this->bounds['ymax']; $y += $cellSize ) {
        $cellQueue->push(new \utility\Cell($x + $h, $y + $h, $h, $polygon));
      }
    }


    // first best guess: polygon centroid
    $this->centroid = $this->getCentroidCell($polygon);

    // second guess: bounding box centroid
    $this->boundsCentroid = new \utility\Cell($this->bounds['xmin'] + $width / 2, $this->bounds['ymin'] + $height / 2, 0, $polygon);


    $bestCell = ( $this->boundsCentroid->d > $this->centroid->d)
      ? $this->boundsCentroid
      : $this->centroid;

    $numProbes = $cellQueue->length();

    while( $cellQueue->length() ) {
      set_time_limit( 60 );
      // pick the most promising cell from the queue
      $cell = $cellQueue->pop();

      // update the best cell if we found a better one
      if( $cell->d > $bestCell->d ) {
        $bestCell = $cell;

        if( $debugCallback ) {
          $debugCallback( sprintf('found best %f after %d probes', round(1e4 * $cell->d) / 1e4, $numProbes) );
        }
      }

      // do not drill down further if there's no chance of a better solution
      if( $cell->max - $bestCell->d <= $precision ) {
        continue;
      }

      // split the cell into four cells
      $h = $cell->h / 2;
      $cellQueue->push(new \utility\Cell($cell->x - $h, $cell->y - $h, $h, $polygon));
      $cellQueue->push(new \utility\Cell($cell->x + $h, $cell->y - $h, $h, $polygon));
      $cellQueue->push(new \utility\Cell($cell->x - $h, $cell->y + $h, $h, $polygon));
      $cellQueue->push(new \utility\Cell($cell->x + $h, $cell->y + $h, $h, $polygon));
      $numProbes += 4;
    }

    if( $debugCallback ) {
      $debugCallback('num probes: ' . $numProbes);
      $debugCallback('best distance: ' . $bestCell->d);
    }

    $this->centerOfMass = [
      'x' => $bestCell->x,
      'y' => $bestCell->y,
      'd' => $bestCell->d
    ];
  }

  // get polygon centroid
  private function getCentroidCell($polygon) {
    $v = &$this->area;

    $x = 0;
    $y = 0;
    $v = 0; // volume

    foreach ( $polygon as $ring ) {
      for( $i = 0, $len = count($ring), $j = $len - 1; $i < $len; $j = $i++ ) {
        $a  = $ring[$i];
        $b  = $ring[$j];
        $f  = $a['x'] * $b['y'] - $b['x'] * $a['y'];
        $x += ($a['x'] + $b['x']) * $f;
        $y += ($a['y'] + $b['y']) * $f;
        $v += $f * 3;
      }
    }

    if( $v === 0 ) {
      return new \utility\Cell($ring[0]['x'], $ring[0]['y'], 0, $polygon);
    }

    return new \utility\Cell($x/$v, $y/$v, 0, $polygon);
  }
}

class Cell {
	public $x;
	public $y;
	public $h;
	public $d;
	public $max;

	public function __construct($x, $y, $h, $polygon) {
		$this->x   = $x; // cell center x
		$this->y   = $y; // cell center y
		$this->h   = $h; // half the cell size
		$this->d   = $this->pointToPolygonDist($x, $y, $polygon);
		$this->max = $this->d + $this->h * M_SQRT2;
	}

  // signed distance from point to polygon outline (negative if point is outside)
  private function pointToPolygonDist( $x, $y, $polygon ) {
    $inside = false;
    $minDistSq = INF;

    for( $k = 0; $k < count($polygon); $k++ ) {
      $ring = $polygon[$k];

      for( $i = 0, $len = count($ring), $j = $len - 1; $i < $len; $j = $i++ ) {
        $a = $ring[$i];
        $b = $ring[$j];

        if(
          ($a['y'] > $y !== $b['y'] > $y) &&
          ($x < ($b['x'] - $a['x']) * ($y - $a['y']) / ($b['y'] - $a['y']) + $a['x'])
        ) {
          $inside = !$inside;
        }

        $minDistSq = min($minDistSq, $this->getSegDistSq($x, $y, $a, $b));
      }

    }

    return $minDistSq === 0 ? 0 : ($inside ? 1 : -1) * sqrt($minDistSq);
  }


  // get squared distance from a point to a segment
  private function getSegDistSq($px, $py, $a, $b) {
    $x  = $a['x'];
    $y  = $a['y'];
    $dx = $b['x'] - $x;
    $dy = $b['y'] - $y;

    if( $dx > 0 || $dy > 0 ) {
      $t = (($px - $x) * $dx + ($py - $y) * $dy) / ($dx * $dx + $dy * $dy);

      if( $t > 1 ) {
        $x = $b['x'];
        $y = $b['y'];
      } else if( $t > 0 ) {
        $x += $dx * $t;
        $y += $dy * $t;
      }
    }

    $dx = $px - $x;
    $dy = $py - $y;

    return $dx * $dx + $dy * $dy;
  }
}

class CellQueue {
	public $queue;

  public function __construct() {
		$this->queue = new \SplPriorityQueue();
	}

	public function push(Cell $cell) {
		$this->queue->insert($cell, $cell->max);
	}

	public function length() {
		return $this->queue->count();
	}

	/** @return Cell */
	public function pop() {
		return $this->queue->extract();
	}
}
