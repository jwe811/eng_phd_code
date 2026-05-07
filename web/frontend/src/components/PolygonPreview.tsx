import React from 'react';
import * as THREE from 'three';
import { api, UofsObject } from '../lib/api';

const colors = [0x2563eb, 0xdc2626, 0x059669, 0x9333ea];

export function PolygonPreview(props: { path: string; index: number }) {
  const mountRef = React.useRef<HTMLDivElement | null>(null);
  const [data, setData] = React.useState<UofsObject | null>(null);
  const [error, setError] = React.useState('');

  React.useEffect(() => {
    let cancelled = false;
    setError('');
    api.uofsObject(props.path, props.index)
      .then((object) => {
        if (!cancelled) setData(object);
      })
      .catch((err) => {
        if (!cancelled) setError(err instanceof Error ? err.message : String(err));
      });
    return () => {
      cancelled = true;
    };
  }, [props.path, props.index]);

  React.useEffect(() => {
    const mount = mountRef.current;
    if (!mount || !data) return;
    mount.innerHTML = '';
    const width = mount.clientWidth || 720;
    const height = mount.clientHeight || 420;
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0xf8fafc);
    const camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(width, height);
    mount.appendChild(renderer.domElement);

    const allPoints = data.polygons.flatMap((poly) => poly.points);
    const center = new THREE.Vector3();
    for (const point of allPoints) center.add(new THREE.Vector3(point[0], point[1], point[2]));
    center.divideScalar(Math.max(1, allPoints.length));
    const maxRadius = Math.max(
      4,
      ...allPoints.map((point) => new THREE.Vector3(point[0], point[1], point[2]).distanceTo(center))
    );

    camera.position.set(center.x + maxRadius * 2.2, center.y + maxRadius * 2.2, center.z + maxRadius * 1.7);
    camera.lookAt(center);
    scene.add(new THREE.AmbientLight(0xffffff, 0.7));
    const light = new THREE.DirectionalLight(0xffffff, 0.7);
    light.position.set(4, 5, 7);
    scene.add(light);

    const grid = new THREE.GridHelper(Math.max(4, maxRadius * 3), Math.max(4, Math.ceil(maxRadius * 3)));
    grid.position.set(center.x, center.y, center.z - 0.02);
    scene.add(grid);

    data.polygons.forEach((poly, index) => {
      const material = new THREE.LineBasicMaterial({ color: colors[index % colors.length], linewidth: 3 });
      const points = poly.points.map((point) => new THREE.Vector3(point[0], point[1], point[2]));
      scene.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
      const nodeMaterial = new THREE.MeshStandardMaterial({ color: colors[index % colors.length] });
      for (const point of points) {
        const node = new THREE.Mesh(new THREE.SphereGeometry(0.055 * maxRadius, 16, 12), nodeMaterial);
        node.position.copy(point);
        scene.add(node);
      }
    });

    let dragging = false;
    let lastX = 0;
    function rotate(delta: number) {
      const offset = camera.position.clone().sub(center);
      offset.applyAxisAngle(new THREE.Vector3(0, 0, 1), delta);
      camera.position.copy(center.clone().add(offset));
      camera.lookAt(center);
    }
    function onPointerDown(event: PointerEvent) {
      dragging = true;
      lastX = event.clientX;
    }
    function onPointerMove(event: PointerEvent) {
      if (!dragging) return;
      rotate((event.clientX - lastX) * 0.01);
      lastX = event.clientX;
      renderer.render(scene, camera);
    }
    function onPointerUp() {
      dragging = false;
    }
    renderer.domElement.addEventListener('pointerdown', onPointerDown);
    renderer.domElement.addEventListener('pointermove', onPointerMove);
    window.addEventListener('pointerup', onPointerUp);
    renderer.render(scene, camera);

    return () => {
      renderer.domElement.removeEventListener('pointerdown', onPointerDown);
      renderer.domElement.removeEventListener('pointermove', onPointerMove);
      window.removeEventListener('pointerup', onPointerUp);
      renderer.dispose();
      mount.innerHTML = '';
    };
  }, [data]);

  return (
    <div>
      {error && <div className="errors">{error}</div>}
      {data?.validation_errors.length ? <div className="errors">{data.validation_errors.join(' ')}</div> : null}
      <div ref={mountRef} className="canvas-wrap" aria-label="3D polygon preview" />
    </div>
  );
}

