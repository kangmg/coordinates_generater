![image](https://github.com/kangmg/coordinates_generater/assets/59556369/4b9979ae-49b0-4d2c-9c27-acd790fe4e5e)

## 사용방법
디렉토리 내에 molecule.sdf란 이름의 파일이 존재해야 합니다.

## 파일 설명
* rvmv.py     : rovolution - inter-fragment move를 수행하고 각 좌표를 coordinates 폴더에 저장함
* dihedral.py : dihedral angle에 대해 rotation을 수행하고 각 xyz좌표를 coordinates 폴더에 저장함

## 생성 파일
* coordinates/*.xyz      : 파일 이름은 dihedral의 경우 회전각.xyz 형태로 저장됩니다. rvmv의 경우 이동거리_moved_공전각_rotated.xyz로 저장됩니다.
* trajectory.xyz         : 파일은 trajectory 형태로 각 파일의 xyz 데이터를 담고 있습니다. Chemcraft와 같은 프로그램을 이용하면 회전 에니메이션 등을 얻을 수 있습니다.
* total_coordinates.xyz  : 모든 xyz파일의 좌표들 담고 있습니다. avogadro와 같은 프로그램을 통해서 실행 결과를 한 눈에 확인할 수 있습니다.
